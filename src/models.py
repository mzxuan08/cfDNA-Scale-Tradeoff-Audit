"""Six fixed-configuration classifiers. DANN uses source domains only."""
from __future__ import annotations

import gc
import math
import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

from .common import seed_everything


def cuda_preflight():
    import torch
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA is required; install a CUDA-enabled PyTorch build and compatible NVIDIA driver')
    return torch


class TorchClassifier(ClassifierMixin, BaseEstimator):
    def __init__(self, domain_adversarial=False, hidden=128, batch_size=128,
                 epochs=80, lr=0.001, seed=2026, require_cuda=True, dropout=.2, weight_decay=.0001):
        self.domain_adversarial = domain_adversarial
        self.hidden, self.batch_size, self.epochs = hidden, batch_size, epochs
        self.lr, self.seed, self.require_cuda = lr, seed, require_cuda
        self.dropout,self.weight_decay=dropout,weight_decay

    def fit(self, X, y, center_train=None):
        import torch
        from torch import nn
        if not 1 <= self.hidden <= 256 or not 1 <= self.batch_size <= 128:
            raise ValueError('Strict GPU constraint: hidden<=256, batch_size<=128')
        if self.epochs < 1:
            raise ValueError('epochs must be positive')
        seed_everything(self.seed)
        torch.set_num_threads(1)
        self.device_ = 'cuda' if self.require_cuda else 'cpu'
        if self.require_cuda:
            cuda_preflight()
            torch.cuda.set_per_process_memory_fraction(0.45)
            torch.cuda.reset_peak_memory_stats()
        x = torch.from_numpy(np.asarray(X, dtype=np.float32))
        labels = torch.from_numpy(np.asarray(y, dtype=np.float32))
        if set(np.unique(y)) != {0, 1}:
            raise ValueError('Both classes required')
        if center_train is None:
            center_train = np.zeros(len(y), dtype=int)
        _, domains = np.unique(center_train, return_inverse=True)
        n_domains = len(np.unique(domains))
        if self.domain_adversarial and n_domains < 2:
            raise ValueError('Source-only DANN requires at least two TRAINING centers')
        domain_labels = torch.from_numpy(domains.astype(np.int64))
        # Keep the full feature matrix in host RAM; send only minibatches to CUDA.
        self.encoder_ = nn.Sequential(nn.Linear(x.shape[1], self.hidden), nn.ReLU(), nn.Dropout(self.dropout),
                                      nn.Linear(self.hidden, self.hidden // 2 or 1), nn.ReLU(), nn.Dropout(self.dropout)).to(self.device_)
        width = self.hidden // 2 or 1
        self.head_ = nn.Linear(width, 1).to(self.device_)
        self.domain_head_ = nn.Linear(width, n_domains).to(self.device_) if self.domain_adversarial else None
        modules = [self.encoder_, self.head_] + ([self.domain_head_] if self.domain_head_ else [])
        parameters = [p for module in modules for p in module.parameters()]
        optimizer = torch.optim.Adam(parameters, lr=self.lr,weight_decay=self.weight_decay)
        generator = torch.Generator().manual_seed(self.seed)
        loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(x, labels, domain_labels),
                                             batch_size=self.batch_size, shuffle=True,
                                             generator=generator, num_workers=0)
        class Reverse(torch.autograd.Function):
            @staticmethod
            def forward(ctx, values, strength):
                ctx.strength = strength
                return values.view_as(values)

            @staticmethod
            def backward(ctx, grad):
                return -ctx.strength * grad, None
        for epoch in range(self.epochs):
            strength = 2 / (1 + math.exp(-10 * epoch / max(self.epochs - 1, 1))) - 1
            for xb, yb, db in loader:
                xb, yb, db = xb.to(self.device_), yb.to(self.device_), db.to(self.device_)
                optimizer.zero_grad(set_to_none=True)
                h = self.encoder_(xb)
                loss = nn.functional.binary_cross_entropy_with_logits(self.head_(h).flatten(), yb)
                if self.domain_adversarial:
                    loss = loss + nn.functional.cross_entropy(self.domain_head_(Reverse.apply(h, strength)), db)
                loss.backward()
                optimizer.step()
                if self.require_cuda and torch.cuda.max_memory_reserved() >= 4 * 1024 ** 3:
                    raise MemoryError('PyTorch allocator exceeded the 4 GiB VRAM budget')
        self.classes_ = np.array([0, 1])
        self.n_features_in_ = x.shape[1]
        self.peak_vram_bytes_ = int(torch.cuda.max_memory_reserved()) if self.require_cuda else 0
        self.encoder_.eval()
        self.head_.eval()
        return self

    def predict_proba(self, X):
        import torch
        x = torch.from_numpy(np.asarray(X, dtype=np.float32))
        values = []
        with torch.no_grad():
            for xb in x.split(self.batch_size):
                values.append(torch.sigmoid(self.head_(self.encoder_(xb.to(self.device_)))).cpu().numpy().ravel())
        p = np.concatenate(values)
        return np.column_stack((1 - p, p)).astype(np.float32)

    def close(self):
        import torch
        for name in ('encoder_', 'head_', 'domain_head_'):
            if hasattr(self, name):
                delattr(self, name)
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


def make_model(name, cfg):
    name={'ElasticNet':'elasticnet','RandomForest':'rf','XGBoost_GPU':'xgboost',
          'LightGBM':'lightgbm','PyTorch_MLP':'mlp','PyTorch_DANN':'dann'}.get(name,name)
    e = cfg['evaluation']
    seed = cfg['seed']
    if name == 'elasticnet':
        return LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5,
                                  C=1.0, max_iter=20000, tol=0.001, random_state=seed)
    if name == 'rf':
        return RandomForestClassifier(n_estimators=300, max_depth=8, min_samples_leaf=2,
                                      n_jobs=e['n_jobs'], random_state=seed)
    if name == 'xgboost':
        if not e['require_cuda']:
            raise ValueError('Production XGBoost requires CUDA by project specification')
        # Independent of torch; validate booster device after fit to catch silent fallback.
        from xgboost import XGBClassifier
        return XGBClassifier(tree_method='hist', device='cuda', n_estimators=300,
                             max_depth=5, learning_rate=0.05, random_state=seed,
                             n_jobs=e['n_jobs'], eval_metric='logloss')
    if name == 'lightgbm':
        from lightgbm import LGBMClassifier
        return LGBMClassifier(n_estimators=200, num_leaves=15, learning_rate=0.05,
                              deterministic=True, force_col_wise=True, verbosity=-1,
                              random_state=seed, n_jobs=e['n_jobs'])
    if name in ('mlp', 'dann'):
        return TorchClassifier(name == 'dann', e['torch_hidden'], e['torch_batch_size'],
                               e['torch_epochs'], e['torch_lr'], seed, e['require_cuda'])
    raise ValueError(f'Unknown model {name}')
