import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from src.evaluation import bh_fdr, delong_test, ece, bootstrap_ci, paired_bootstrap_delta, metric_values
from src.variance_audit import partial_permanova


def test_delong_against_hand_computed_placements():
    y = [1, 1, 0, 0]
    r = delong_test(y, [.8, .4, .6, .2], [.9, .8, .2, .1])
    assert r['auc_a'] == .75
    assert r['auc_b'] == 1
    assert r['variance'] == pytest.approx(.125)
    assert r['p'] == pytest.approx(2 * norm.sf(.25 / np.sqrt(.125)))


def test_delong_ties_and_identical_predictions():
    p = [.5, .5, .5, .5]
    r = delong_test([0, 1, 0, 1], p, p)
    assert r['auc_a'] == .5 and r['p'] == 1 and r['delta'] == 0


def test_bh_known_values_and_missingness():
    np.testing.assert_allclose(bh_fdr([.01, .04, .03, .002, np.nan]), [.02, .04, .04, .008, np.nan], equal_nan=True)
    with pytest.raises(ValueError):
        bh_fdr([1.1])


def test_calibration_boundaries():
    assert ece([0, 1], [0, 1]) == 0
    assert ece([1, 0], [0, 1]) == 1
    assert np.isnan(metric_values([1, 1], [.2, .3])['auroc'])


def test_bootstrap_is_reproducible_and_paired():
    f = pd.DataFrame(dict(patient_id=[f'p{i}' for i in range(60)], center=np.repeat(['A', 'B', 'C'], 20),
                          label=[0, 1] * 30, probability=[.2, .8] * 30))
    a, b = bootstrap_ci(f, repeats=60), bootstrap_ci(f, repeats=60)
    assert a == b and a['low'] == a['high'] == 1
    f['random_probability'] = f.probability
    f['loco_probability'] = f.probability
    r = paired_bootstrap_delta(f, 60)
    assert r['delta_random_minus_loco'] == r['ci_low'] == r['ci_high'] == 0


def test_permanova_strong_center_signal_is_detected():
    rng = np.random.default_rng(2026)
    centers = np.repeat([0, 1, 2], 20)
    y = np.tile([0, 1], 30)
    X = (centers[:, None] * 10 + rng.normal(size=(60, 3))).astype('float32')
    reduced = np.column_stack([np.ones(60), y])
    center_design = pd.get_dummies(centers, drop_first=True).to_numpy(dtype=float)
    r = partial_permanova(X, reduced, center_design, y, permutations=49)
    assert r['partial_R2'] > .9 and r['p'] <= .04


def test_permanova_confounded_design_rejected():
    y = np.array([0, 0, 1, 1])
    reduced = np.column_stack([np.ones(4), y])
    with pytest.raises(ValueError, match='not identifiable'):
        partial_permanova(np.arange(8).reshape(4, 2), reduced, y[:, None], y)
