import numpy as np
import pandas as pd
import pytest

from src.debiasing import TrainMedianScaler, InFoldComBat, FoldPreprocessor
from src.evaluation import make_splits, assert_split


@pytest.mark.parametrize('n_centers', [3, 4, 6])
def test_loco_every_fold_purges_cross_center_patient_aliases(n_centers):
    frame = pd.DataFrame([dict(sample_id=f's{c}_{i}', patient_id=f'p{c}_{i}', center=f'C{c}', label=i % 2)
                          for c in range(n_centers) for i in range(12)])
    # The same control patient appears in every center; all other identities differ.
    frame.loc[frame.index % 12 == 0, 'patient_id'] = 'shared-control'
    folds = make_splits(frame, 'loco')
    assert len(folds) == n_centers
    tested = []
    for train, test in folds:
        train_centers, test_centers = frame.iloc[train].center, frame.iloc[test].center
        train_patients, test_patients = frame.iloc[train].patient_id, frame.iloc[test].patient_id
        assert set(train_centers).intersection(set(test_centers)) == set()
        assert set(train_patients).intersection(set(test_patients)) == set()
        assert len(set(test_centers)) == 1
        tested.extend(test)
    assert sorted(tested) == list(range(len(frame)))


@pytest.mark.parametrize('method', ['raw', 'lsn', 'icc', 'lsn_icc', 'combat'])
def test_transform_never_mutates_locked_mean_variance_or_feature_mask(method):
    import pickle
    from sklearn.base import BaseEstimator, TransformerMixin
    # Complete cross-center measurements make ICC estimable in this test only.
    base = np.arange(1, 9, dtype='float32')
    values = np.column_stack([base, base ** 2, base ** 3, base + .3])
    X = np.repeat(values, 2, axis=0)
    centers = np.tile(['A', 'B'], 8)
    patients = np.repeat([f'p{i}' for i in range(8)], 2)
    names = ['ratio__chr1:0-5', 'ratio__chr1:5-10', 'ratio__chr1:10-15', 'length__80']
    transformer = FoldPreprocessor(method, names, threshold=.1, min_subjects=5).fit(
        X, np.repeat([0, 1], 8), centers, patients)
    assert isinstance(transformer, (BaseEstimator, TransformerMixin))
    snapshot = pickle.dumps(transformer)
    mean, variance = transformer.scaler_.mean_.copy(), transformer.scaler_.var_.copy()
    support = transformer.support_.copy()
    test = np.array([[1e4, np.nan, 2e4, 1e5], [2e4, 1e4, 4e4, 2e5]], dtype='float32')
    output = transformer.transform(test, center_test=['A', 'B']) if method == 'combat' else transformer.transform(test)
    assert output.dtype == np.float32 and np.isfinite(output).all()
    np.testing.assert_array_equal(transformer.scaler_.mean_, mean)
    np.testing.assert_array_equal(transformer.scaler_.var_, variance)
    np.testing.assert_array_equal(transformer.support_, support)
    assert pickle.dumps(transformer) == snapshot


@pytest.mark.parametrize('mode', ['random', 'loco'])
def test_patient_and_center_isolation(cohort, mode):
    cohort.loc[12, 'patient_id'] = cohort.loc[0, 'patient_id']
    for train, test in make_splits(cohort, mode, 3):
        a, b = cohort.iloc[train], cohort.iloc[test]
        assert not set(a.patient_id) & set(b.patient_id)
        if mode == 'loco':
            assert not set(a.center) & set(b.center)


def test_split_guard_catches_patient_leak(cohort):
    cohort.loc[12, 'patient_id'] = cohort.loc[0, 'patient_id']
    with pytest.raises(ValueError, match='Patient leakage'):
        assert_split(cohort, [0, 1], [12, 13])


def test_scaler_uses_train_median_and_keeps_empty_columns():
    model = TrainMedianScaler().fit([[1, np.nan], [3, np.nan]])
    before = model.median_.copy()
    result = model.transform([[np.nan, np.nan], [1000, 10]])
    np.testing.assert_equal(before, model.median_)
    np.testing.assert_allclose(result, [[0, 0], [998, 10]])
    assert result.dtype == np.float32


def test_combat_rejects_unseen_center_without_fitting():
    rng = np.random.default_rng(2026)
    train = rng.normal(size=(20, 8)).astype('float32')
    train[10:] += 3
    f = InFoldComBat().fit(train, center_train=['A'] * 10 + ['B'] * 10)
    before = {c: tuple(x.copy() for x in p) for c, p in f.params_.items()}
    with pytest.raises(ValueError, match='unseen center'):
        f.transform(train[:1], center_test=['C'])
    for c, params in before.items():
        for old, current in zip(params, f.params_[c]):
            np.testing.assert_array_equal(old, current)
    out = f.transform(train, center_test=['A'] * 10 + ['B'] * 10)
    assert np.isfinite(out).all()
    assert np.abs(out[10:].mean() - out[:10].mean()) < np.abs(train[10:].mean() - train[:10].mean())


def test_fold_preprocessing_invariant_to_other_test_samples():
    cols = ['ratio__chr1:0-5', 'ratio__chr1:5-10', 'length__80']
    f = FoldPreprocessor('lsn', cols).fit([[1, 3, .2], [2, 4, .4]])
    alone = f.transform([[10, 20, .3]])
    mixed = f.transform([[10, 20, .3], [100000, 0, .9]])
    np.testing.assert_array_equal(alone[0], mixed[0])


def test_loco_rejects_single_class_training(cohort):
    cohort.loc[cohort.center.ne('A'), 'label'] = 0
    with pytest.raises(ValueError, match='Training fold'):
        make_splits(cohort, 'loco', 3)
