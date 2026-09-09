"""Time-series estimators of irreversibility (entropy-production proxies) plus surrogate nulls.

Every estimator here returns 0 in expectation for a time-reversible process, so the sign and
magnitude are only interpretable against surrogates that destroy irreversibility while preserving
everything else. Nothing in this module should be used without the surrogate comparison.

Estimators
  ehlers(x)        skewness of first differences. Zero if the increment distribution is symmetric.
  trev(x, lag)     third-order time-reversibility statistic <(x_{t+L}-x_t)^3>, normalised.
                   NOTE at lag=1 this is algebraically identical to ehlers (both reduce to the
                   skewness of first differences), so the default lag is 4 to probe a genuinely
                   different timescale and keep the two estimators independent.
  perm_kl(x, m)    KL divergence between ordinal-pattern distributions of the forward and
                   time-reversed series. Closest in spirit to a Roldan-Parrondo path-KL estimate
                   of entropy production, and the least sensitive to the marginal distribution.

Surrogates
  shuffle_surrogate       destroys all temporal structure   -> every estimator must go to 0
  phase_surrogate         preserves the power spectrum, destroys phase relations (and hence
                          nonlinear/irreversible structure) -> every estimator must go to 0

CRITICAL: all of these are sample-size biased. perm_kl especially, because it estimates a KL
divergence over m! bins from a finite sample. Never compare windows of different length -- match
sample COUNTS. See `zscore_vs_surrogates`, which returns a bias-corrected score.
"""
from __future__ import annotations
import numpy as np
import math
from itertools import permutations


# ----------------------------------------------------------------- estimators
def ehlers(x):
    """Skewness of first differences; 0 for a time-reversible process."""
    d = np.diff(np.asarray(x, float))
    s = d.std()
    return 0.0 if s == 0 else float(np.mean(d ** 3) / s ** 3)


def trev(x, lag=4):
    """Third-order time-reversibility statistic, normalised to be scale-free."""
    x = np.asarray(x, float)
    d = x[lag:] - x[:-lag]
    v = np.mean(d ** 2)
    return 0.0 if v == 0 else float(np.mean(d ** 3) / v ** 1.5)


_PERM_INDEX = {}


def _perm_codes(x, m):
    """Ordinal (permutation) pattern codes of embedding order m."""
    x = np.asarray(x, float)
    n = len(x) - m + 1
    if n <= 0:
        return np.empty(0, int)
    if m not in _PERM_INDEX:
        _PERM_INDEX[m] = {p: i for i, p in enumerate(permutations(range(m)))}
    idx = _PERM_INDEX[m]
    win = np.lib.stride_tricks.sliding_window_view(x, m)
    order = np.argsort(win, axis=1, kind="stable")
    return np.fromiter((idx[tuple(r)] for r in order), int, count=n)


def perm_kl(x, m=3):
    """KL(forward || reversed) over ordinal patterns. Zero iff the pattern statistics are
    time-symmetric. Laplace-smoothed so the divergence is always finite."""
    f = _perm_codes(x, m)
    r = _perm_codes(np.asarray(x, float)[::-1], m)
    if len(f) < 2 or len(r) < 2:
        return 0.0
    k = math.factorial(m)
    pf = (np.bincount(f, minlength=k) + 0.5); pf = pf / pf.sum()
    pr = (np.bincount(r, minlength=k) + 0.5); pr = pr / pr.sum()
    return float(np.sum(pf * np.log(pf / pr)))


def epr_markov(x, nbins=5, lag=1):
    """Plug-in entropy-production rate from the empirical transition matrix.

    Symbolise the series into `nbins` equiprobable states, count transitions, and evaluate the
    antisymmetric flux contribution

        EPR = (1/2) * sum_ij (N_ij - N_ji)/N * ln(N_ij / N_ji)

    Unlike the skewness- and KL-based statistics, this is proportional to the true entropy
    production rate, which is what a scaling-exponent test requires. It is strongly
    sample-size biased (it needs nbins^2 transition counts), so nbins must stay small and the
    surrogate correction is mandatory.
    """
    x = np.asarray(x, float)
    if len(x) < 4 * nbins:
        return 0.0
    edges = np.quantile(x, np.linspace(0, 1, nbins + 1)[1:-1])
    s = np.searchsorted(edges, x)
    a, b = s[:-lag], s[lag:]
    N = np.zeros((nbins, nbins))
    np.add.at(N, (a, b), 1.0)
    tot = N.sum()
    if tot == 0:
        return 0.0
    out = 0.0
    for i in range(nbins):
        for j in range(i + 1, nbins):
            f, r = N[i, j], N[j, i]
            if f > 0 and r > 0:
                out += (f - r) / tot * np.log(f / r)
    return float(out)


ESTIMATORS = {"ehlers": ehlers, "trev4": lambda x: trev(x, 4),
              "perm_kl": lambda x: perm_kl(x, 3), "epr_markov": lambda x: epr_markov(x, 5)}


# ----------------------------------------------------------------- surrogates
def shuffle_surrogate(x, rng):
    return rng.permutation(np.asarray(x, float))


def phase_surrogate(x, rng):
    """Randomise Fourier phases, preserving the power spectrum (and hence the autocorrelation)."""
    x = np.asarray(x, float)
    n = len(x)
    F = np.fft.rfft(x - x.mean())
    ph = rng.uniform(0, 2 * np.pi, len(F))
    ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    return np.fft.irfft(np.abs(F) * np.exp(1j * ph), n=n) + x.mean()


SURROGATES = {"shuffle": shuffle_surrogate, "phase": phase_surrogate}


def zscore_vs_surrogates(x, estimator, surrogate="phase", nsurr=50, seed=0):
    """Bias-corrected irreversibility: (value - surrogate mean) / surrogate sd.

    The surrogate mean absorbs the finite-sample bias, which is why raw estimator values
    must never be compared across windows of different length.
    """
    rng = np.random.default_rng(seed)
    f = ESTIMATORS[estimator] if isinstance(estimator, str) else estimator
    g = SURROGATES[surrogate]
    obs = f(x)
    null = np.array([f(g(x, rng)) for _ in range(nsurr)])
    sd = null.std(ddof=1)
    return dict(value=obs, null_mean=float(null.mean()), null_sd=float(sd),
                z=float((obs - null.mean()) / sd) if sd > 0 else 0.0,
                excess=float(obs - null.mean()))
