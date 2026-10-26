"""Reference-derived annotations and exploratory, source-traceable mechanism tests."""
import gzip
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,fisher_exact
from .common import write_json,write_parquet,sha256
from .reference import load_reference_windows
from .stream_extractor import remote_operation,HashReader,stream_gc_reference,make_bins
from .cohort import evidence_download
from .study import load_clinical,export_table
from .debiasing import ArmLocalLogRatio
from .evaluation import bh_fdr


def interval_fractions(chunks,bins):
    gathered={c:[] for c in bins.chrom.unique()}
    for frame in chunks:
        for chrom,group in frame[frame.chrom.isin(gathered)].groupby('chrom',sort=False):
            gathered[chrom].append(group[['start','end']].to_numpy(dtype=np.int64))
    values=np.zeros(len(bins),dtype=np.float64)
    for chrom,parts in gathered.items():
        if not parts:continue
        a=np.concatenate(parts);a=a[np.argsort(a[:,0],kind='stable')]
        running=np.maximum.accumulate(a[:,1]);starts=np.r_[0,np.flatnonzero(a[1:,0]>running[:-1])+1]
        begin=a[starts,0];end=np.maximum.reduceat(a[:,1],starts);prefix=np.r_[0,np.cumsum(end-begin)]
        def integral(x):
            index=np.searchsorted(end,x,side='right');safe=np.minimum(index,len(begin)-1)
            return prefix[index]+np.where(index<len(begin),np.maximum(0,x-begin[safe]),0)
        group=bins[bins.chrom.eq(chrom)]
        values[group.index]=(integral(group.end.to_numpy())-integral(group.start.to_numpy()))/(group.end-group.start)
    return values.astype('float32')


def stream_interval_annotation(url,bins,cfg,columns):
    def consume(response):
        reader=HashReader(response.raw)
        with gzip.GzipFile(fileobj=reader) as handle:
            frames=pd.read_csv(handle,sep='\t',header=None,usecols=columns,chunksize=200000,comment='#')
            def normalize():
                for frame in frames:
                    frame.columns=['chrom','start','end']
                    yield frame
            fractions=interval_fractions(normalize(),bins)
        return fractions,{'url':url,'compressed_sha256':reader.digest.hexdigest(),'compressed_bytes':reader.count}
    return remote_operation(url,cfg['network'],consume)


def build_annotations(cfg):
    ref=Path(cfg['paths']['reference']);bins=load_reference_windows(cfg).reset_index(drop=True)
    path=ref/'hg19_window_annotations.parquet'
    if path.exists() and path.with_suffix('.json').exists():
        stamp=json.loads(path.with_suffix('.json').read_text())
        if stamp.get('compartment_parser_version')==2 and stamp['sha256']==sha256(path) and stamp['window_sha256']==sha256(ref/cfg['reference']['window_table']):return pd.read_parquet(path)
    result=bins.copy();provenance={}
    gcpath=ref/'hg19_5mb_gc.parquet'
    if gcpath.exists():gc=pd.read_parquet(gcpath)
    else:
        full=make_bins(ref/'hg19.chrom.sizes',cfg['features']['bin_size'],cfg['features']['chromosomes'])
        gc,source=stream_gc_reference(cfg['reference']['gc_urls']['hg19'],full,cfg)
        write_parquet(gc,gcpath);write_json(source,gcpath.with_suffix('.json'))
    result=result.merge(gc[['bin_id','gc_fraction','called_fraction']],on='bin_id',validate='one_to_one')
    provenance['gc']=json.loads(gcpath.with_suffix('.json').read_text())
    print('GC annotation complete',flush=True)
    for key,column,columns in [('repeatmasker','repeat_fraction',[5,6,7]),('segmental_duplications','segdup_fraction',[1,2,3])]:
        cache=ref/(key+'_fractions.parquet')
        if cache.exists():data=pd.read_parquet(cache);fraction=data[column].to_numpy();source=json.loads(cache.with_suffix('.json').read_text())
        else:
            fraction,source=stream_interval_annotation(cfg['study']['annotations'][key],result,cfg,columns)
            write_parquet(pd.DataFrame({'bin_id':result.bin_id,column:fraction}),cache);write_json(source,cache.with_suffix('.json'))
        result[column]=fraction;provenance[key]=source
        print(key+' annotation complete',flush=True)
    # Compartment labels are a GM12878 reference, not a tumor-specific assay.
    try:
        url=cfg['study']['annotations']['compartments'];file=evidence_download(url,ref/'GM12878_subcompartments.bed.gz',cfg['network'])
        frame=pd.read_csv(file,sep='\t',header=None,comment='#',keep_default_na=False)
        if frame.shape[1]<4:raise ValueError('Subcompartment table missing labels')
        frame=frame.iloc[:,:4];frame.columns=['chrom','start','end','compartment']
        known=frame.compartment.astype(str).str.match(r'^[AB][0-9]*$')
        if not known.any():raise ValueError('No explicit A/B subcompartment labels')
        aframe=frame[known & frame.compartment.str.startswith('A')][['chrom','start','end']]
        result['A_compartment_fraction']=interval_fractions([aframe],result)
        coverage=interval_fractions([frame.loc[known,['chrom','start','end']]],result)
        result['compartment_coverage_fraction']=coverage
        result.loc[coverage==0,'A_compartment_fraction']=np.nan
        provenance['compartments']={'url':url,'sha256':sha256(file),'unclassified_intervals':int((~known).sum()),
                  'reference':'GM12878 lymphoblastoid, hg19; not tumor-specific; A fraction is covered bp / full window, not an imputed B label'}
    except Exception as exc:
        result['A_compartment_fraction']=np.float32(np.nan)
        provenance['compartments']={'unavailable':str(exc),'no_compartment_values_invented':True}
    write_parquet(result,path)
    write_json({'sha256':sha256(path),'window_sha256':sha256(ref/cfg['reference']['window_table']),
                'compartment_parser_version':2,'sources':provenance},path.with_suffix('.json'))
    return result


def biology_analysis(cfg):
    started=time.perf_counter()
    import statsmodels.api as sm
    from .plotting import clinical_figure4
    annotations=build_annotations(cfg)
    meta,R,_,bins=load_clinical(cfg,True)
    L=ArmLocalLogRatio(bins,cfg['study']['radius']).fit(R).transform(R)
    destination=Path(cfg['paths']['results']);folder=destination/'tables'
    partition=pd.read_csv(folder/'Bin_Variance_Partition.csv').query('family=="raw_ratio"')
    data=annotations.merge(partition[['bin_id','eta2_batch','eta2_bio','vulnerability_ratio']],on='bin_id',validate='one_to_one')
    data['gc_extremeness']=np.abs(data.gc_fraction-data.gc_fraction.median())
    data['log_telomere_distance']=np.log1p(data.distance_to_telomere_bp)
    data['log_centromere_distance']=np.log1p(data.distance_to_centromere_bp)
    correlations=[]
    variables=['gc_fraction','gc_extremeness','log_telomere_distance','log_centromere_distance','repeat_fraction','segdup_fraction','A_compartment_fraction']
    for variable in variables:
        valid=data[[variable,'eta2_batch']].dropna()
        if len(valid)<5:continue
        r,p=spearmanr(valid[variable],valid.eta2_batch)
        correlations.append(dict(annotation=variable,rho=float(r),p=float(p),n=len(valid),interpretation='exploratory_genomic_association_not_causal_proof'))
    correlation=pd.DataFrame(correlations);correlation['q']=bh_fdr(correlation.p)
    export_table(correlation,cfg,'Genomic_Spearman_Associations')
    predictors=['gc_extremeness','log_telomere_distance','log_centromere_distance','repeat_fraction','segdup_fraction']
    design=data[predictors].astype(float);design=(design-design.mean())/design.std()
    model=sm.OLS(data.eta2_batch,sm.add_constant(design)).fit(cov_type='cluster',cov_kwds={'groups':data.chrom})
    ci=model.conf_int();regression=pd.DataFrame({'term':model.params.index,'beta':model.params.to_numpy(),
        'p':model.pvalues.to_numpy(),'ci_low':ci.iloc[:,0].to_numpy(),'ci_high':ci.iloc[:,1].to_numpy()})
    export_table(regression,cfg,'Genomic_Multivariable_Regression')
    effect=[]
    for center in sorted(meta.center.unique()):
        h=meta.center.eq(center)&meta.label.eq(0);c=meta.center.eq(center)&meta.label.eq(1)
        if min(h.sum(),c.sum())<2:continue
        for name,x in [('Raw_Ratio',R),('LSN',L)]:
            sd=np.sqrt(((h.sum()-1)*np.var(x[h],axis=0,ddof=1)+(c.sum()-1)*np.var(x[c],axis=0,ddof=1))/(h.sum()+c.sum()-2))
            d=(np.mean(x[c],axis=0)-np.mean(x[h],axis=0))/np.maximum(sd,1e-8)
            effect.extend(dict(bin_id=bin_id,center=center,pipeline=name,cohen_d=float(value),healthy_n=int(h.sum()),cancer_n=int(c.sum())) for bin_id,value in zip(bins.bin_id,d))
    effects=pd.DataFrame(effect);export_table(effects,cfg,'Disease_Effect_Sizes')
    selected=[]
    for file in sorted((destination/'clinical_cv_checkpoints').glob('selection_loco_*_Proposed_LSN_ICC.parquet')):
        selected.append(pd.read_parquet(file))
    if len(selected)!=3:raise ValueError('Three completed LOCO masks required before biological filtering interpretation')
    selection=pd.concat(selected).groupby('bin_id',as_index=False).agg(retained_fraction=('retained','mean'),healthy_center_icc=('healthy_center_icc','mean'))
    data=data.merge(selection,on='bin_id',validate='one_to_one');data['removed_majority']=data.retained_fraction<.5
    flags={'GC_extreme_deciles':data.gc_fraction.lt(data.gc_fraction.quantile(.1))|data.gc_fraction.gt(data.gc_fraction.quantile(.9)),
           'Within_10Mb_telomere':data.distance_to_telomere_bp.lt(10000000),
           'Within_10Mb_centromere':data.distance_to_centromere_bp.lt(10000000),
           'High_repeat_fraction':data.repeat_fraction.gt(data.repeat_fraction.median()),
           'High_segmental_duplication':data.segdup_fraction.gt(data.segdup_fraction.median()),
           'Hypothesized_CNA_arms_1q_7p_8q':data.chromosome_arm.isin(['chr1q','chr7p','chr8q'])}
    if data.A_compartment_fraction.notna().all():flags['GM12878_A_compartment_majority']=data.A_compartment_fraction.gt(.5)
    enrichment=[]
    for name,flag in flags.items():
        removed=data.removed_majority
        contingency=np.array([[np.sum(removed&flag),np.sum(removed&~flag)],[np.sum(~removed&flag),np.sum(~removed&~flag)]])
        _,p=fisher_exact(contingency)
        adjusted=contingency.astype(float)+(.5 if (contingency==0).any() else 0)
        a,b,c,d=adjusted.ravel();odds=a*d/(b*c);se=np.sqrt(np.sum(1/adjusted))
        enrichment.append(dict(region=name,odds_ratio=odds,ci_low=np.exp(np.log(odds)-1.96*se),
                     ci_high=np.exp(np.log(odds)+1.96*se),p=float(p),a=int(contingency[0,0]),b=int(contingency[0,1]),
                     c=int(contingency[1,0]),d=int(contingency[1,1]),zero_correction=(contingency==0).any()))
    enrichment=pd.DataFrame(enrichment);enrichment['q']=bh_fdr(enrichment.p)
    export_table(enrichment,cfg,'Filtered_Window_Enrichment')
    importance=pd.read_csv(folder/'Clinical_Feature_Importance.csv')
    importance=importance[(importance['mode']=='loco')&(importance.pipeline=='Proposed_LSN_ICC')&(importance.model=='elasticnet')]
    ranks=importance.pivot(index='bin_id',columns='fold',values='rank')
    concordance=[]
    for a in range(3):
        for b in range(a+1,3):
            pair=ranks[[a,b]].dropna();r,p=spearmanr(pair[a],pair[b]);concordance.append(dict(fold_a=a,fold_b=b,rho=float(r),p=float(p),shared_windows=len(pair)))
    export_table(pd.DataFrame(concordance),cfg,'Feature_Rank_Concordance')
    export_table(data,cfg,'Genomic_Interpretability_Windows')
    clinical_figure4(data,correlation,enrichment,meta,R,L,bins,importance,cfg)
    paired=effects.pivot(index=['bin_id','center'],columns='pipeline',values='cohen_d').dropna()
    summary={'annotation_windows':len(data),'removed_majority_windows':int(data.removed_majority.sum()),
             'effect_pairs':len(paired),'fraction_effect_magnitude_larger_after_lsn':float((paired.LSN.abs()>paired.Raw_Ratio.abs()).mean()),
             'median_absolute_d_raw':float(paired.Raw_Ratio.abs().median()),'median_absolute_d_lsn':float(paired.LSN.abs().median()),
             'regression_R2':float(model.rsquared),'regression_SE':'chromosome_cluster_robust_22_clusters',
             'interpretation':'Reference associations and local contrasts do not establish tumor-specific nucleosome depletion or causal sequence mechanisms',
             'missing_compartment_annotation':bool(data.A_compartment_fraction.isna().any()),
             'one_control_study_excluded_from_Cohens_d':True,
             'runtime_seconds':round(time.perf_counter()-started,3),
             'independent_counts_per_study':meta.groupby('center').agg(patients=('patient_id','size'),cancer=('label','sum')).to_dict('index'),
             'annotation_missing_fraction':data.isna().mean().to_dict(),
             'clinical_metadata_missing_fraction':meta.isna().mean().to_dict()}
    write_json(summary,destination/'reports/step_06_biology_summary.json')
    (destination/'reports/step_06_biology_report.md').write_text('# Genomic interpretability\n\n'+correlation.to_markdown(index=False)+'\n\n'+enrichment.to_markdown(index=False)+'\n\n'+json.dumps(summary,indent=2),encoding='utf-8')
    return summary
