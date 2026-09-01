"""A4 control: circulation of a SINGLE colloidal cluster, matched to the single-axoneme analysis.

Referee objection: the colloid circulation null was measured on system-averaged network
observables over ~10^4 particles, while the cilia measurement uses the two-mode shape space of a
SINGLE axoneme. In an isotropic, statistically homogeneous suspension most global scalar pairs
have vanishing signed area by symmetry, and phase-incoherent circulation across many clusters
averages away -- precisely the incoherence documented in Sec. 3.5 (C sqrt(N) ~ const). The colloid
may therefore fail for reasons of ensemble averaging and observable choice rather than because its
irreversibility dies under coarse-graining.

This applies the identical estimator at a matched level of description: one tracked cluster,
represented by a low-order shape descriptor, PCA-reduced to two modes, signed area rate against
phase-randomised surrogates. If a single active cluster also gives |z| < 1.5, the claim survives
and is much stronger; if it does not, the conclusion changes.

Shape descriptor: the deviatoric gyration tensor of the cluster (Gxx-Gyy, 2Gxy) plus the
normalised third and fourth mass multipoles, which is the closest available analogue of the
axoneme tangent-angle representation for a compact 2-D object -- translation- and
rotation-referenced, and low-order.
"""
from __future__ import annotations
import glob
import numpy as np
import networkx as nx
from analyze import contact_graph


def cluster_shape(pos, members, L):
    """Low-order shape descriptor of one cluster, translation-referenced."""
    p = pos[members]
    d = p - p[0]
    d -= L * np.round(d / L)              # unwrap about a reference member
    d -= d.mean(0)                        # centre of mass frame
    x, y = d[:, 0], d[:, 1]
    r2 = (x ** 2 + y ** 2).mean()
    if r2 <= 0:
        return None
    Gxx, Gyy, Gxy = (x * x).mean(), (y * y).mean(), (x * y).mean()
    z = (x + 1j * y) / np.sqrt(r2)
    m3, m4 = np.mean(z ** 3), np.mean(z ** 4)
    return np.array([(Gxx - Gyy) / r2, 2 * Gxy / r2,
                     m3.real, m3.imag, m4.real, m4.imag])


def area_rate(x, y, dt):
    x = (x - x.mean()) / (x.std() + 1e-30)
    y = (y - y.mean()) / (y.std() + 1e-30)
    return 0.5 * np.sum(x[:-1] * y[1:] - y[:-1] * x[1:]) / (len(x) * dt)


def phase_surr(x, rng):
    F = np.fft.rfft(x - x.mean()); ph = rng.uniform(0, 2 * np.pi, len(F)); ph[0] = 0
    if len(x) % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(F) * np.exp(1j * ph), n=len(x)) + x.mean()


def track_and_measure(f, min_size=30, nsurr=60):
    z = np.load(f)
    sn, sv, L = z["snaps"], z["svec"], float(z["Lbox"])
    dts = float(z["dt_snap"])
    lo = len(sn) // 2
    # track the largest cluster by membership overlap across consecutive snapshots
    prev, series = None, []
    for t in range(lo, len(sn)):
        G = contact_graph(sn[t], sv, L)
        comps = [c for c in nx.connected_components(G) if len(c) >= min_size]
        if not comps:
            break
        if prev is None:
            cur = max(comps, key=len)
        else:
            ov = [len(c & prev) for c in comps]
            if max(ov) < 0.5 * len(prev):
                break                      # identity lost: stop rather than switch objects
            cur = comps[int(np.argmax(ov))]
        s = cluster_shape(sn[t], np.fromiter(cur, int), L)
        if s is None:
            break
        series.append(s); prev = cur
    if len(series) < 40:
        return None
    S = np.array(series)
    S = S - S.mean(0)
    u, sv_, _ = np.linalg.svd(S, full_matrices=False)
    a1, a2 = u[:, 0] * sv_[0], u[:, 1] * sv_[1]
    obs = area_rate(a1, a2, dts)
    rng = np.random.default_rng(0)
    null = np.array([area_rate(phase_surr(a1, rng), phase_surr(a2, rng), dts)
                     for _ in range(nsurr)])
    sd = null.std(ddof=1)
    return dict(n=len(series), size=len(prev), area=obs, null=null.mean(),
                z=(obs - null.mean()) / sd if sd > 0 else 0.0)


if __name__ == "__main__":
    print("=== single-cluster circulation, matched to the single-axoneme analysis ===")
    print("  cilia (single axoneme, 2 shape modes): median |z| = 3.2, 92% above |z| = 2\n")
    print("   condition                      n_clusters   median |z|   frac |z|>2   tracked frames")
    for tag, pat in [("chi=0   (reciprocal)", "data_paper/*_x0_*N4000*.npz"),
                     ("chi=1   (nonreciprocal)", "data_paper/*_x1_*N4000*seed[123].npz"),
                     ("chi=1.5 (nonreciprocal)", "data_paper/*_x1.5*N4000*cont1.npz")]:
        R = [r for r in (track_and_measure(f) for f in sorted(glob.glob(pat))) if r]
        if not R:
            print(f"   {tag:28s}  -- no trackable cluster --"); continue
        zz = np.abs([r["z"] for r in R])
        print(f"   {tag:28s} {len(R):8d}     {np.median(zz):8.2f}     {np.mean(zz > 2):7.2f}"
              f"       {int(np.median([r['n'] for r in R])):6d}")
