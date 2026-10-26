import gzip
import importlib.util
import io
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import sys
import threading

import numpy as np
import pandas as pd
import pytest
import yaml

from src.common import ROOT, load_config, write_parquet, load_dataset
from src.cohort import select_cohort, stratified_cap
from src.stream_extractor import aggregate_fragment_batches, aggregate_fragments, histogram_quantiles, Blacklist


def test_histogram_quartiles_equal_expanded_numpy_linear():
    lengths = np.array([80, 80, 100, 150, 151, 220, 400])
    hist = np.bincount(lengths - 80, minlength=321)[None, :]
    np.testing.assert_allclose(histogram_quantiles(hist)[0], np.quantile(lengths, [.25, .5, .75]))
    assert np.isnan(histogram_quantiles(np.zeros((1, 321), dtype=int))).all()


def test_vectorized_parser_does_not_shift_into_filtered_window_gaps(extract_params):
    params = {**extract_params, 'chromosomes': ['chr1'], 'normalize_chr_prefix': True}
    bins = pd.DataFrame([dict(chrom='chr1', start=0, end=5000000, chrom_size_bp=15000000, bin_id='chr1:0-5000000'),
                         dict(chrom='chr1', start=10000000, end=15000000, chrom_size_bp=15000000, bin_id='chr1:10000000-15000000')])
    source = io.BytesIO(b'1\t0\t100\t60\t+\n1\t0\t150\t60\t+\n1\t0\t151\t60\t+\n'
                         b'1\t5000000\t5000220\t60\t+\n1\t10000000\t10000400\t60\t+\n'
                         b'1\t0\t150\t29\t+\nchrX\t0\t150\t60\t+\n')
    values, qc = aggregate_fragment_batches(source, bins, params, Blacklist(), chunk_rows=2)
    assert qc['accepted'] == 5 and qc['outside_retained_windows'] == 1
    assert qc['low_mapq'] == 1 and qc['excluded_chrom'] == 1
    assert values['short_count__chr1:0-5000000'] == 2
    assert values['long_count__chr1:0-5000000'] == 1
    assert values['total_count__chr1:0-5000000'] == 3
    assert values['iqr__chr1:0-5000000'] == pytest.approx(np.subtract(*np.quantile([100,150,151],[.75,.25])))
    assert values['total_count__chr1:10000000-15000000'] == 0
    assert values['length__220'] == pytest.approx(.2)
    assert sum(v for k,v in values.items() if k.startswith('length__')) == pytest.approx(1)
    assert sum(v for k,v in values.items() if k.startswith('length10__')) == pytest.approx(1)


def test_vectorized_and_scalar_extraction_agree_with_realistic_filters(reference_bins, extract_params):
    rng=np.random.default_rng(2026)
    params={**extract_params,'chromosomes':['chr1'],'normalize_chr_prefix':True}
    starts=rng.integers(0,9_999_000,size=500)
    lengths=rng.integers(70,420,size=500)
    mapqs=rng.choice([0,29,30,60],size=500)
    lines=[f'chr1\t{s}\t{s+l}\t{q}\t+\n' for s,l,q in zip(starts,lengths,mapqs)]
    lines+=['chr1\t10000\t10150\t60\t+\n','chr1\t10100\t10251\t60\t+\n',
            'chr1\t10200\t10351\t60\t+\n']
    mask=Blacklist();mask.intervals['chr1']=np.array([[10000,10200],[4999990,5000010]])
    expected,expected_qc=aggregate_fragments(lines,reference_bins,params,mask)
    actual,qc=aggregate_fragment_batches(io.BytesIO(''.join(lines).encode()),reference_bins,params,mask,chunk_rows=17)
    for name,value in expected.items():
        np.testing.assert_allclose(actual[name],value,rtol=1e-6,equal_nan=True)
    assert qc['accepted']==expected_qc['accepted']


def raw_fixture():
    return pd.DataFrame([dict(type='Healthy' if i%2==0 else 'Breast', patient_id=f'patient_{c}_{i}',
          specimen_id=f'specimen_{c}_{i}', timepoint=1, bam_id=f'EE{1000+c*100+i}',study=f'Study {c}',
          sample_type='plasma', stage='I' if i%2 else None,seq_depth=1.,gender='female')
          for c in range(3) for i in range(10)])


def test_stratified_selection_never_relabels_nonhealthy_controls():
    cfg = yaml.safe_load((ROOT/'config/config.yaml').read_text())
    cfg['cohort'].update(min_samples=20,max_samples=30,min_per_center=8)
    frame = raw_fixture()
    benign = frame.iloc[[0]].copy()
    benign['type'],benign['bam_id'],benign['patient_id'],benign['specimen_id']='Cirrhosis','EE9999','benign','benign'
    selected,screening = select_cohort(pd.concat([frame,benign],ignore_index=True),cfg)
    assert len(selected)==30 and screening.loc[screening.bam_id.eq('EE9999'),'reason'].iloc[0]=='not_cancer_or_healthy'
    one=stratified_cap(selected,24,2026);two=stratified_cap(selected,24,2026)
    assert one.sample_id.tolist()==two.sample_id.tolist()
    assert one.groupby('center').label.nunique().eq(2).all()


def test_parallel_stream_checkpoint_resume_qc_and_publication(tmp_path, monkeypatch):
    cfg = yaml.safe_load((ROOT/'config/config.yaml').read_text())
    cfg['cohort'].update(min_samples=20,max_samples=30,min_per_center=8)
    cfg['features']['chromosomes']=['chr1']
    cfg['extraction'].update(chunk_rows=2,min_effective_fragments=1)
    (tmp_path/'config').mkdir()
    config=tmp_path/'config/config.yaml'
    class Handler(BaseHTTPRequestHandler):
        calls=0
        def do_GET(self):
            Handler.calls+=1
            # One genuine QC failure in the TEST-only fixture; no fake production input.
            payload=(b'chr1\t0\t150\t29\t+\n' if 'EE1000' in self.path else
                     b'chr1\t0\t150\t60\t+\nchr1\t0\t151\t60\t-\n')
            self.send_response(200);self.end_headers();self.wfile.write(gzip.compress(payload))
        def log_message(self,*args):pass
    server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    cfg['metadata']['public_s3_base']=f'http://127.0.0.1:{server.server_port}'
    config.write_text(yaml.safe_dump(cfg))
    loaded=load_config(config)
    ref=Path(loaded['paths']['reference'])
    (ref/'hg19.chrom.sizes').write_text('chr1\t5000000\n')
    (ref/'hg19.blacklist.bed.gz').write_bytes(gzip.compress(b''))
    bins=pd.DataFrame([dict(chrom='chr1',start=0,end=5000000,chrom_size_bp=5000000,bin_id='chr1:0-5000000')])
    write_parquet(bins,ref/cfg['reference']['window_table'])
    sys.path.insert(0,str(ROOT/'scripts'))
    spec=importlib.util.spec_from_file_location('step02_test',ROOT/'scripts/02_fetch_and_extract.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    monkeypatch.setattr(module,'published_metadata',lambda _:raw_fixture())
    monkeypatch.setattr(module,'load_reference_windows',lambda _:bins)
    monkeypatch.setattr(sys,'argv',[str(ROOT/'scripts/02_fetch_and_extract.py'),'--config',str(config),'--skip-live-probe'])
    try:
        # A partial run has truly unknown statistics in the unprocessed studies.
        monkeypatch.setattr(sys,'argv',[str(ROOT/'scripts/02_fetch_and_extract.py'),'--config',str(config),
                                        '--skip-live-probe','--max-new-samples','6'])
        module.main()
        import json
        progress=json.loads((Path(loaded['paths']['results'])/'reports/step_02_progress.json').read_text())
        assert progress['completed']==6 and progress['pending']==24
        assert any(r['median_effective_fragments'] is None for r in progress['study_counts'])
        monkeypatch.setattr(sys,'argv',[str(ROOT/'scripts/02_fetch_and_extract.py'),'--config',str(config),'--skip-live-probe'])
        module.main()
        assert Handler.calls==30
        final=Path(loaded['paths']['processed'])/cfg['extraction']['final_file']
        assert len(pd.read_parquet(final))==29
        table=pd.read_csv(Path(loaded['paths']['results'])/'tables/Table1_Cohort_Characteristics.csv')
        assert table.qc_exclusions.sum()==1 and table.pending_samples.sum()==0
        module.main()
        assert Handler.calls==30,'Completed sample streams must not be downloaded again'
        meta,X=load_dataset(load_config(config))
        assert len(meta)==29 and X.to_numpy().dtype==np.float32
    finally:
        server.shutdown();server.server_close();thread.join()
