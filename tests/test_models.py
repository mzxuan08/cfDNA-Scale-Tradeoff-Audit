import numpy as np
import pytest

from src.models import make_model, TorchClassifier


@pytest.mark.parametrize('name', ['elasticnet', 'rf', 'lightgbm', 'xgboost', 'mlp', 'dann'])
def test_six_models_fit_and_predict_on_test_only_data(name):
    if name in ('mlp', 'dann'):
        torch = pytest.importorskip('torch')
        if not torch.cuda.is_available():
            pytest.skip('CUDA not available; production CUDA requirement is not relaxed')
    if name in ('xgboost', 'lightgbm'):
        pytest.importorskip(name)
    rng = np.random.default_rng(2026)
    y = np.tile([0, 1], 30)
    X = rng.normal(size=(60, 8)).astype('float32')
    X[:, 0] += y * 2
    cfg = dict(seed=2026, evaluation=dict(n_jobs=1, require_cuda=True, torch_hidden=16,
               torch_batch_size=16, torch_epochs=2, torch_lr=.001))
    model = make_model(name, cfg)
    if name in ('rf', 'xgboost', 'lightgbm'):
        model.set_params(n_estimators=5)
    try:
        if isinstance(model, TorchClassifier):
            model.fit(X[:48], y[:48], center_train=np.repeat(['A', 'B'], 24))
        else:
            model.fit(X[:48], y[:48])
        p = model.predict_proba(X[48:])
        assert p.shape == (12, 2) and np.isfinite(p).all()
        np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-6)
        if name == 'xgboost':
            import json
            device = json.loads(model.get_booster().save_config())['learner']['generic_param']['device']
            assert device.startswith('cuda'), 'CUDA fallback is not a successful GPU test'
        if isinstance(model, TorchClassifier):
            assert model.peak_vram_bytes_ < 4 * 1024 ** 3
    finally:
        if isinstance(model, TorchClassifier):
            model.close()
