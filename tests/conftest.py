from pathlib import Path
import sys
import pytest
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


@pytest.fixture
def cohort():
    return pd.DataFrame([dict(sample_id=f's{c}{i}', patient_id=f'p{c}{i}', center=c, label=i % 2)
                         for c in ['A', 'B', 'C'] for i in range(12)])


@pytest.fixture
def reference_bins():
    return pd.DataFrame([('chr1', 0, 5000000, 'chr1:0-5000000'),
                         ('chr1', 5000000, 10000000, 'chr1:5000000-10000000')],
                        columns=['chrom', 'start', 'end', 'bin_id'])


@pytest.fixture
def extract_params():
    return dict(histogram=[80, 400], short=[100, 150], long=[151, 220], bin_size=5000000,
                pseudocount=.5, min_long_count=1, min_accepted_fragments=1, mapq_min=30)
