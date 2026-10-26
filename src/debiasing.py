"""Train-only transformations with explicit unseen-center and ICC design contracts."""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted


def array(X):
    a = np.asarray(X, dtype=np.float32)
    if a.ndim != 2 or np.isinf(a).any():
        raise ValueError('Expected finite-or-NaN 2D float32 feature matrix')
    return a


class TrainMedianScaler(BaseEstimator, TransformerMixin):
    """Preserve all-missing columns as zero; means/scales learned on training only."""
    def fit(self, X_train, y_train=None, center_train=None):
        X = array(X_train)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            median = np.nanmedian(X, axis=0)
        self.median_ = np.nan_to_num(median, nan=0).astype(np.float32)
        clean = np.where(np.isnan(X), self.median_, X)
        self.mean_ = clean.mean(axis=0, dtype=np.float64).astype(np.float32)
        self.var_ = clean.var(axis=0, dtype=np.float64).astype(np.float32)
        self.scale_ = clean.std(axis=0, dtype=np.float64).astype(np.float32)
        self.scale_[self.scale_ < 1e-7] = 1
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X_test):
        check_is_fitted(self, 'median_')
        X = array(X_test)
        if X.shape[1] != self.n_features_in_:
            raise ValueError('Feature count mismatch')
        return ((np.where(np.isnan(X), self.median_, X) - self.mean_) / self.scale_).astype(np.float32)


class LocalSelfNormalizer(BaseEstimator, TransformerMixin):
    """R_ij / median(R_ik) for finite neighbors on the same chromosome.

    Window includes focal bin and at most radius adjacent bins on each side.
    Histogram columns pass through. No test population statistics are used.
    NaN focal values and zero/undefined local medians remain NaN.
    """
    def __init__(self, feature_names, radius=2):
        self.feature_names = feature_names
        self.radius = radius

    def fit(self, X_train, y_train=None, center_train=None):
        X = array(X_train)
        if X.shape[1] != len(self.feature_names) or self.radius < 1:
            raise ValueError('Invalid LSN feature names/radius')
        groups = {}
        for i, name in enumerate(self.feature_names):
            if name.startswith('ratio__'):
                chrom, interval = name.split('__', 1)[1].split(':')
                start = int(interval.split('-')[0])
                groups.setdefault(chrom, []).append((start, i))
        self.neighbors_ = {}
        for values in groups.values():
            idx = [i for _, i in sorted(values)]
            for k, j in enumerate(idx):
                self.neighbors_[j] = idx[max(0, k - self.radius):k + self.radius + 1]
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X_test):
        check_is_fitted(self, 'neighbors_')
        X = array(X_test)
        if X.shape[1] != self.n_features_in_:
            raise ValueError('Feature count mismatch')
        out = X.copy()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            for j, indices in self.neighbors_.items():
                denominator = np.nanmedian(X[:, indices], axis=1)
                valid = np.isfinite(denominator) & (denominator > 0)
                out[:, j] = np.nan
                out[valid, j] = X[valid, j] / denominator[valid]
        return out


def icc_absolute_agreement(panel):
    """ICC(A,1)/ICC(2,1), complete subjects x random measurement centers.

    (MS_subject - MS_error) /
    (MS_subject + (k-1)*MS_error + k*(MS_center-MS_error)/n).
    Negative ICCs are retained. Undefined constant features return NaN.
    """
    x = np.asarray(panel, dtype=np.float64)
    if x.ndim != 2 or min(x.shape) < 2 or not np.isfinite(x).all():
        raise ValueError('ICC requires >=2 complete subjects and >=2 centers')
    n, k = x.shape
    mean = x.mean()
    subject, center = x.mean(axis=1), x.mean(axis=0)
    ms_subject = k * np.sum((subject - mean) ** 2) / (n - 1)
    ms_center = n * np.sum((center - mean) ** 2) / (k - 1)
    residual = x - subject[:, None] - center[None, :] + mean
    ms_error = np.sum(residual ** 2) / ((n - 1) * (k - 1))
    denominator = ms_subject + (k - 1) * ms_error + k * (ms_center - ms_error) / n
    return float((ms_subject - ms_error) / denominator) if denominator > 1e-15 else float('nan')


class ICCDesignError(ValueError):
    pass


class ICCStabilityFilter(BaseEstimator, TransformerMixin):
    def __init__(self, threshold=0.75, min_subjects=5):
        self.threshold = threshold
        self.min_subjects = min_subjects

    def fit(self, X_train, y_train=None, center_train=None, subject_train=None):
        X = array(X_train)
        if center_train is None or subject_train is None:
            raise ICCDesignError('ICC requires explicit repeated-subject and center identifiers')
        centers, subjects = np.asarray(center_train), np.asarray(subject_train)
        if len(centers) != len(X) or len(subjects) != len(X):
            raise ICCDesignError('ICC metadata length mismatch')
        frame = pd.DataFrame({'subject': subjects, 'center': centers, 'row': np.arange(len(X))})
        if frame.duplicated(['subject', 'center']).any():
            raise ICCDesignError('Choose one measurement per subject-center cell before ICC')
        lookup = frame.pivot(index='subject', columns='center', values='row').dropna()
        if lookup.shape[1] < 2 or len(lookup) < self.min_subjects:
            raise ICCDesignError(f'ICC requires >=2 training centers and >= {self.min_subjects} '
                                 f'complete shared subjects; observed {lookup.shape}')
        rows = lookup.to_numpy(dtype=int)
        scores = []
        for j in range(X.shape[1]):
            panel = X[rows, j]
            panel = panel[np.isfinite(panel).all(axis=1)]
            scores.append(icc_absolute_agreement(panel) if len(panel) >= self.min_subjects else np.nan)
        self.icc_ = np.asarray(scores)
        self.support_ = np.isfinite(self.icc_) & (self.icc_ >= self.threshold)
        self.subjects_ = lookup.index.to_list()
        if not self.support_.any():
            raise ICCDesignError('No features meet the prespecified ICC threshold; no fallback selection')
        return self

    def transform(self, X_test):
        check_is_fitted(self, 'support_')
        X = array(X_test)
        if X.shape[1] != len(self.support_):
            raise ValueError('Feature count mismatch')
        return X[:, self.support_]


class InFoldComBat(BaseEstimator, TransformerMixin):
    """Parametric empirical-Bayes ComBat, without outcome covariates.

    Features borrow strength for location/scale priors. Train-only parameters.
    Transform needs centers (or an index level named center). An unseen center
    raises: estimating its parameters from test samples would change the protocol.
    """
    def __init__(self, max_iter=1000, tolerance=1e-5):
        self.max_iter = max_iter
        self.tolerance = tolerance

    def fit(self, X_train, y_train=None, center_train=None):
        X = np.asarray(X_train, dtype=np.float64)
        if center_train is None or not np.isfinite(X).all():
            raise ValueError('ComBat needs center IDs and train-only imputed finite input')
        centers = np.asarray(center_train)
        unique, counts = np.unique(centers, return_counts=True)
        if len(unique) < 2 or min(counts) < 2 or X.shape[1] < 2:
            raise ValueError('ComBat requires >=2 centers, >=2 samples/center, >=2 features')
        means = np.stack([X[centers == c].mean(axis=0) for c in unique])
        self.mean_ = np.average(means, weights=counts, axis=0)
        residual = X.copy()
        for c, mean in zip(unique, means):
            residual[centers == c] -= mean
        self.scale_ = np.sqrt(np.maximum(np.mean(residual ** 2, axis=0), 1e-8))
        Z = (X - self.mean_) / self.scale_
        self.params_ = {}
        for c in unique:
            z = Z[centers == c]
            n = len(z)
            gamma_hat = z.mean(axis=0)
            delta_hat = np.maximum(z.var(axis=0, ddof=1), 1e-8)
            gamma_bar = gamma_hat.mean()
            t2 = max(gamma_hat.var(ddof=1), 1e-8)
            m, v = delta_hat.mean(), max(delta_hat.var(ddof=1), 1e-8)
            a, b = (2 * v + m * m) / v, (m * v + m ** 3) / v
            gamma, delta = gamma_hat.copy(), delta_hat.copy()
            for _ in range(self.max_iter):
                g = (t2 * n * gamma_hat + delta * gamma_bar) / (t2 * n + delta)
                d = (0.5 * np.sum((z - g) ** 2, axis=0) + b) / (n / 2 + a - 1)
                change = max(np.max(np.abs(g - gamma) / (np.abs(gamma) + 1e-8)),
                             np.max(np.abs(d - delta) / (delta + 1e-8)))
                gamma, delta = g, np.maximum(d, 1e-8)
                if change < self.tolerance:
                    break
            else:
                raise RuntimeError('ComBat posterior iteration did not converge')
            self.params_[c] = (gamma, delta)
        return self

    def transform(self, X_test, center_test=None):
        check_is_fitted(self, 'params_')
        if center_test is None and isinstance(X_test, pd.DataFrame) and 'center' in X_test.index.names:
            center_test = X_test.index.get_level_values('center').to_numpy()
        if center_test is None:
            raise ValueError('ComBat transform requires center_test or center-indexed DataFrame')
        X = np.asarray(X_test, dtype=np.float64)
        centers = np.asarray(center_test)
        if len(centers) != len(X) or not np.isfinite(X).all():
            raise ValueError('Invalid ComBat test data')
        if set(centers) - self.params_.keys():
            raise ValueError('ComBat cannot transform an unseen center without test-center fitting')
        Z = (X - self.mean_) / self.scale_
        for c in np.unique(centers):
            gamma, delta = self.params_[c]
            Z[centers == c] = (Z[centers == c] - gamma) / np.sqrt(delta)
        return (Z * self.scale_ + self.mean_).astype(np.float32)


class FoldPreprocessor(BaseEstimator, TransformerMixin):
    def __init__(self, method, feature_names, radius=2, threshold=0.75, min_subjects=5):
        self.method, self.feature_names = method, feature_names
        self.radius, self.threshold, self.min_subjects = radius, threshold, min_subjects

    def fit(self, X_train, y_train=None, center_train=None, subject_train=None):
        if self.method not in ('raw', 'lsn', 'icc', 'lsn_icc', 'combat'):
            raise ValueError(f'Unknown method {self.method}')
        X = array(X_train)
        self.lsn_ = self.icc_filter_ = self.combat_ = None
        if self.method in ('lsn', 'lsn_icc'):
            self.lsn_ = LocalSelfNormalizer(self.feature_names, self.radius).fit(X)
            X = self.lsn_.transform(X)
        self.support_ = np.ones(X.shape[1], dtype=bool)
        if self.method in ('icc', 'lsn_icc'):
            self.icc_filter_ = ICCStabilityFilter(self.threshold, self.min_subjects).fit(
                X, y_train, center_train, subject_train)
            self.support_ = self.icc_filter_.support_
            X = self.icc_filter_.transform(X)
        self.scaler_ = TrainMedianScaler().fit(X, y_train, center_train)
        X = self.scaler_.transform(X)
        if self.method == 'combat':
            self.combat_ = InFoldComBat().fit(X, center_train=center_train)
        return self

    def transform(self, X_test, center_test=None):
        check_is_fitted(self, 'scaler_')
        X = array(X_test)
        if self.lsn_:
            X = self.lsn_.transform(X)
        if self.icc_filter_:
            X = self.icc_filter_.transform(X)
        X = self.scaler_.transform(X)
        if self.combat_:
            X = self.combat_.transform(X, center_test)
        return X


class ArmLocalLogRatio(BaseEstimator, TransformerMixin):
    """Sample-wise log2 contrast to up to k OTHER windows on each side, same arm.

    epsilon is relative to the sample's genomic median. This preserves global
    multiplicative scale equivariance; additive and arbitrary nonlinear invariance
    are neither mathematically implied nor claimed.
    """
    def __init__(self, bins, radius=3, epsilon=1e-6):
        self.bins, self.radius, self.epsilon = bins, radius, epsilon

    def fit(self, X_train, y_train=None, center_train=None):
        X = array(X_train)
        if self.radius < 1 or X.shape[1] != len(self.bins):
            raise ValueError('Invalid LSN radius or genomic schema')
        self.neighbors_ = {}
        for _, group in self.bins.reset_index(drop=True).groupby(['chrom','arm'], sort=False):
            ids = group.sort_values('start').index.to_list()
            for position, idx in enumerate(ids):
                self.neighbors_[idx] = ids[max(0,position-self.radius):position] + ids[position+1:position+1+self.radius]
        return self

    def transform(self, X_test):
        check_is_fitted(self, 'neighbors_')
        X = array(X_test)
        if X.shape[1] != len(self.bins) or np.nanmin(X) < 0:
            raise ValueError('LSN requires nonnegative ratios and the original bin schema')
        eps = self.epsilon * np.maximum(np.nanmedian(X, axis=1), np.finfo(np.float32).tiny)
        output = np.full_like(X, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            for idx, neighbors in self.neighbors_.items():
                if neighbors:
                    baseline = np.nanmedian(X[:,neighbors], axis=1)
                    output[:,idx] = np.log2((X[:,idx]+eps)/(baseline+eps))
        return output.astype('float32')


def healthy_center_components(X, y, centers):
    """Unequal-group one-way random center ICC(1) among different healthy people.

    This measures center clustering, not repeat-measurement reliability.
    n0=(N-sum(n_j^2)/N)/(K-1); variance_center=max((MSB-MSW)/n0,0).
    """
    x = np.asarray(X, dtype=np.float64)[np.asarray(y)==0]
    center = np.asarray(centers)[np.asarray(y)==0]
    groups = [x[center==c] for c in np.unique(center)]
    if len(groups)<2 or len(x)<=len(groups):
        raise ValueError('Healthy center variance needs >=2 centers and positive residual df')
    n = np.array([len(g) for g in groups], dtype=float)
    means = np.stack([g.mean(axis=0) for g in groups])
    grand = x.mean(axis=0)
    ssb = np.sum(n[:,None]*(means-grand)**2, axis=0)
    ssw = sum(np.sum((g-g.mean(axis=0))**2, axis=0) for g in groups)
    msb, msw = ssb/(len(groups)-1), ssw/(len(x)-len(groups))
    n0 = (len(x)-np.sum(n*n)/len(x))/(len(groups)-1)
    between = np.maximum((msb-msw)/n0,0)
    icc = np.divide(between, between+msw, out=np.zeros_like(between), where=(between+msw)>0)
    return icc, between, msw, n.astype(int)


class HealthyCenterStabilityFilter(BaseEstimator, TransformerMixin):
    def __init__(self, exclusion_fraction=.25):
        self.exclusion_fraction = exclusion_fraction

    def fit(self, X_train, y_train=None, center_train=None):
        x = array(X_train)
        if y_train is None or center_train is None or not 0<=self.exclusion_fraction<=.5:
            raise ValueError('Training labels/centers and alpha in [0,.5] required')
        self.imputer_ = TrainMedianScaler().fit(x)
        clean = np.where(np.isnan(x),self.imputer_.median_,x)
        self.icc_, self.between_, self.within_, self.healthy_counts_ = healthy_center_components(clean,y_train,center_train)
        y = np.asarray(y_train)
        signal = np.abs(clean[y==1].mean(axis=0)-clean[y==0].mean(axis=0))
        self.snr_ = signal/(np.sqrt(np.maximum(self.within_,0))+1e-8)
        # High healthy center clustering first; low SNR breaks equal-ICC ties.
        order = np.lexsort((np.arange(x.shape[1]), self.snr_, -self.icc_))
        self.support_ = np.ones(x.shape[1],dtype=bool)
        self.support_[order[:int(np.ceil(x.shape[1]*self.exclusion_fraction))]] = False
        return self

    def transform(self, X_test):
        check_is_fitted(self,'support_')
        return array(X_test)[:,self.support_]


class ClinicalPreprocessor(BaseEstimator, TransformerMixin):
    """The four requested clinical pipelines plus prespecified ablations.

    InFold_ComBat is explicitly transductive for unseen centers: a robust location/
    scale mapping uses the unlabeled target batch, without modifying fitted state.
    The other three pipelines never use other test subjects' distributions.
    """
    def __init__(self, pipeline, bins, radius=3, exclusion_fraction=.25, epsilon=1e-6):
        self.pipeline,self.bins,self.radius = pipeline,bins,radius
        self.exclusion_fraction,self.epsilon = exclusion_fraction,epsilon

    def _sample_transform(self,X):
        x=array(X)
        if self.pipeline=='Hybrid_MultiScale':
            global_z=(x-x.mean(axis=1,keepdims=True))/np.maximum(x.std(axis=1,keepdims=True),1e-8)
            return np.concatenate([global_z,self.lsn_.transform(x)],axis=1).astype('float32')
        if self.lsn_ is not None:
            return self.lsn_.transform(x)
        if self.pipeline=='Global_Zscore':
            return ((x-x.mean(axis=1,keepdims=True))/np.maximum(x.std(axis=1,keepdims=True),1e-8)).astype('float32')
        return x

    def fit(self, X_train, y_train=None, center_train=None):
        from sklearn.preprocessing import RobustScaler, StandardScaler
        allowed=['Raw_Ratio','Global_Zscore','InFold_ComBat','Proposed_LSN_ICC','LSN_only','ICC_only','Hybrid_MultiScale']
        if self.pipeline not in allowed:
            raise ValueError(f'Unknown clinical pipeline {self.pipeline}')
        self.lsn_=None
        if self.pipeline in ('Proposed_LSN_ICC','LSN_only','Hybrid_MultiScale'):
            self.lsn_=ArmLocalLogRatio(self.bins,self.radius,self.epsilon).fit(X_train)
        x=self._sample_transform(X_train)
        self.filter_=None
        self.support_=np.ones(x.shape[1],dtype=bool)
        if self.pipeline in ('Proposed_LSN_ICC','ICC_only'):
            self.filter_=HealthyCenterStabilityFilter(self.exclusion_fraction).fit(x,y_train,center_train)
            self.support_=self.filter_.support_
            x=self.filter_.transform(x)
        self.imputer_=TrainMedianScaler().fit(x)
        clean=np.where(np.isnan(x),self.imputer_.median_,x)
        self.combat_=None
        if self.pipeline=='InFold_ComBat':
            self.combat_=InFoldComBat().fit(clean,center_train=center_train)
            transformed=self.combat_.transform(clean,center_test=center_train)
            self.source_median_=np.median(transformed,axis=0)
            self.source_iqr_=np.maximum(np.subtract(*np.percentile(transformed,[75,25],axis=0)),1e-6)
            clean=transformed
        self.scaler_=(RobustScaler() if self.pipeline in ('Proposed_LSN_ICC','LSN_only','ICC_only') else StandardScaler()).fit(clean)
        self.train_centers_=np.unique(center_train) if center_train is not None else np.array([])
        return self

    def transform(self, X_test, center_test=None):
        check_is_fitted(self,'scaler_')
        if center_test is None and hasattr(X_test,'index') and 'center' in X_test.index.names:
            center_test=X_test.index.get_level_values('center').to_numpy()
        x=self._sample_transform(X_test)
        if self.filter_ is not None:x=self.filter_.transform(x)
        x=np.where(np.isnan(x),self.imputer_.median_,x)
        if self.combat_ is not None:
            if center_test is None:raise ValueError('ComBat needs explicit center metadata')
            centers=np.asarray(center_test)
            mapped=np.empty_like(x)
            for c in np.unique(centers):
                mask=centers==c
                if c in self.combat_.params_:
                    mapped[mask]=self.combat_.transform(x[mask],center_test=centers[mask])
                else:
                    # Allowed user-requested target-batch adaptation; no target labels.
                    median=np.median(x[mask],axis=0)
                    iqr=np.maximum(np.subtract(*np.percentile(x[mask],[75,25],axis=0)),1e-6)
                    mapped[mask]=(x[mask]-median)/iqr*self.source_iqr_+self.source_median_
            x=mapped
        return self.scaler_.transform(x).astype('float32')

    def fit_transform(self,X,y=None,**kwargs):
        return self.fit(X,y,**kwargs).transform(X,center_test=kwargs.get('center_train'))
