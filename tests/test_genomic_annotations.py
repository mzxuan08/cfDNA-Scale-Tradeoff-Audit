import numpy as np
import pandas as pd
from src.interpretability import interval_fractions


def test_annotation_union_does_not_double_count_nested_or_cross_chunk_intervals():
    bins=pd.DataFrame({'chrom':['chr1','chr1','chr2','chr3'],
                       'start':[0,10,0,0],'end':[10,20,10,10]})
    first=pd.DataFrame({'chrom':['chr1','chr1','chr2'],'start':[2,4,0],'end':[8,6,3]})
    second=pd.DataFrame({'chrom':['chr1','chr1','chr2','chrX'],
                         'start':[6,18,3,0],'end':[13,25,10,100]})
    # chr1 union [2,13) U [18,25); adjacent chr2 intervals cover the full bin.
    np.testing.assert_allclose(interval_fractions([first,second],bins),[.8,.5,1.,0.])
