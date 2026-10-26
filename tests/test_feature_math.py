import numpy as np
import pytest

from src.debiasing import LocalSelfNormalizer, icc_absolute_agreement, ICCStabilityFilter, ICCDesignError
from src.stream_extractor import aggregate_fragments, Blacklist


def test_lsn_chromosome_boundaries_and_sample_independence():
    names = ['ratio__chr1:0-5', 'ratio__chr1:5-10', 'ratio__chr2:0-5', 'length__80']
    x = np.array([[1, 3, 100, .2], [2, 6, 200, .3]], dtype='float32')
    model = LocalSelfNormalizer(names, 1).fit(x[:1])
    out = model.transform(x)
    np.testing.assert_allclose(out[:, :3], [[.5, 1.5, 1], [.5, 1.5, 1]])
    np.testing.assert_equal(out[:, 3], x[:, 3])
    np.testing.assert_allclose(out[1], model.transform(x[1:])[0])
    assert out.dtype == np.float32


def test_lsn_undefined_zero_and_nan():
    names = ['ratio__chr1:0-5', 'ratio__chr1:5-10']
    model = LocalSelfNormalizer(names, 1).fit([[0, 0]])
    assert np.isnan(model.transform([[0, 0]])).all()
    np.testing.assert_allclose(model.transform([[np.nan, 2]]), [[np.nan, 1]], equal_nan=True)


def test_icc_absolute_agreement_penalizes_systematic_center_shift():
    assert icc_absolute_agreement([[1, 2], [2, 3], [3, 4]]) == pytest.approx(2 / 3)
    assert icc_absolute_agreement([[1, 1], [2, 2], [3, 3]]) == pytest.approx(1)
    assert np.isnan(icc_absolute_agreement([[1, 1], [1, 1]]))
    assert icc_absolute_agreement([[1, 3], [2, 2], [3, 1]]) < 0


def test_icc_requires_repeated_subjects_not_diagnosis_groups():
    with pytest.raises(ICCDesignError, match='shared subjects'):
        ICCStabilityFilter(min_subjects=2).fit([[1], [2], [3], [4]], [0, 0, 1, 1],
                                               ['A', 'A', 'B', 'B'], ['p1', 'p2', 'p3', 'p4'])


def test_icc_selects_true_reliable_feature():
    x = [[1, 1], [1, 6], [2, 1], [2, 6], [3, 1], [3, 6]]
    f = ICCStabilityFilter(.75, 3).fit(x, center_train=['A', 'B'] * 3, subject_train=['p1', 'p1', 'p2', 'p2', 'p3', 'p3'])
    np.testing.assert_equal(f.support_, [True, False])
    assert f.transform([[4, 100]]).tolist() == [[4]]


def test_extraction_inclusive_bounds_histogram_and_midpoint(reference_bins, extract_params):
    lengths = [80, 99, 100, 150, 151, 220, 221, 400, 401]
    lines = [f'chr1\t0\t{length}\t60\t+\n' for length in lengths]
    lines += ['chr1\t4999950\t5000100\t60\t+\n', 'chr1\t5000000\t5000151\t60\t-\n']
    values, qc = aggregate_fragments(lines, reference_bins, extract_params, Blacklist())
    assert qc['accepted'] == 10
    assert qc['length_outside'] == 1
    assert values['short_count__chr1:0-5000000'] == 2
    assert values['long_count__chr1:0-5000000'] == 2
    assert values['ratio__chr1:0-5000000'] == 1
    assert values['short_count__chr1:5000000-10000000'] == 1
    assert values['long_count__chr1:5000000-10000000'] == 1
    assert sum(v for k, v in values.items() if k.startswith('length__')) == pytest.approx(1, abs=1e-6)
    assert values['length__400'] == pytest.approx(.1)


def test_blacklist_overlap_is_half_open_and_merged(tmp_path):
    path = tmp_path / 'blacklist.bed'
    path.write_text('chr1\t100\t200\nchr1\t150\t250\n')
    b = Blacklist(path)
    assert b.overlaps('chr1', 249, 251)
    assert not b.overlaps('chr1', 250, 300)
    assert not b.overlaps('chr1', 50, 100)


def test_zero_long_is_nan_not_artificial_signal(reference_bins, extract_params):
    values, _ = aggregate_fragments(['chr1\t0\t150\t60\t+'], reference_bins, extract_params, Blacklist())
    assert np.isnan(values['ratio__chr1:0-5000000'])


def test_explicit_chromosome_prefix_normalization(reference_bins, extract_params):
    extract_params['normalize_chr_prefix'] = True
    values, qc = aggregate_fragments(['1\t0\t150\t60\t+', '1\t0\t151\t60\t-'], reference_bins, extract_params, Blacklist())
    assert qc['renamed_chromosome'] == qc['accepted'] == 2
    assert values['ratio__chr1:0-5000000'] == 1


def test_wrong_genome_or_malformed_stream_fails(reference_bins, extract_params):
    with pytest.raises(ValueError, match='reference length'):
        aggregate_fragments(['chr1\t0\t20000000\t60\t+'], reference_bins, extract_params, Blacklist())
    with pytest.raises(ValueError, match='expected'):
        aggregate_fragments(['chr1\t0\t100'], reference_bins, extract_params, Blacklist())
