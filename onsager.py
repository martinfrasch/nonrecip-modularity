"""Onsager reciprocity test: two independent nonreciprocal channels.

A uniform field on one species does NOT work as a second drive: at zero field the system is
isotropic, so the species drift vanishes for every chi and the cross-coefficient is zero by
symmetry rather than by physics. The test would read 0 == 0.

Instead we open a SECOND nonreciprocal channel, in the collision force, with its own parameter
chi2. Its antisymmetric part has the same structure as the EHD one -- equal and same-signed on
both members of a pair, i.e. a pair co-propulsion:

    F_col on i:   k(ssum-r) rhat  +  chi2 * h_ij * k(ssum-r) rhat
    F_col on j:  -k(ssum-r) rhat  +  chi2 * h_ij * k(ssum-r) rhat      h_ij = (s_i-s_j)/(s_i+s_j)

Entropy production is then bilinear in two internal thermodynamic forces,

    T * EPR = chi * J1  +  chi2 * J2 ,
    J1 = sum_i f1_i . v_i    (f1 = EHD nonreciprocal force at unit chi)
    J2 = sum_i f2_i . v_i    (f2 = collision nonreciprocal force at unit chi2)

so chi, chi2 are the forces and J1, J2 their conjugate fluxes. Onsager reciprocity is then

    L_12 = dJ1/dchi2   ==   dJ2/dchi = L_21

Both channels drive the same internal degrees of freedom, so the cross-coefficients are
generically nonzero and the test has content.

Central differences are taken with COMMON RANDOM NUMBERS (identical seed, identical starting
configuration), so each derivative is a paired difference and the discretisation bias that
plagued the raw entropy-production estimator largely cancels.
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp, os
import numpy as np
from numba import njit
from simulate import SIGMA, RCUT


@njit(cache=True, fastmath=True)
def _run2(pos, svec, l2, l4, nsteps, dt, Lbox, alpha, chi, chi2, sigma, seed):
    np.random.seed(seed)
    N = pos.shape[0]
    ncell = int(Lbox / RCUT); csize = Lbox / ncell
    head = np.full(ncell * ncell, -1, np.int64); nxt = np.full(N, -1, np.int64)
    F = np.zeros((N, 2)); F1 = np.zeros((N, 2)); F2 = np.zeros((N, 2))
    half = 0.5 * Lbox
    nstd = np.empty(N)
    for i in range(N):
        nstd[i] = sigma * np.sqrt(dt / svec[i])
    accJ1 = 0.0; accJ2 = 0.0
    for step in range(nsteps):
        head[:] = -1
        for i in range(N):
            c = (int(pos[i, 0] / csize) % ncell) * ncell + (int(pos[i, 1] / csize) % ncell)
            nxt[i] = head[c]; head[c] = i
        F[:] = 0.0; F1[:] = 0.0; F2[:] = 0.0
        for cx in range(ncell):
            for cy in range(ncell):
                i = head[cx * ncell + cy]
                while i >= 0:
                    for dx in range(-1, 2):
                        for dy in range(-1, 2):
                            jc = ((cx + dx) % ncell) * ncell + ((cy + dy) % ncell)
                            j = head[jc]
                            while j >= 0:
                                if j > i:
                                    rx = pos[i, 0] - pos[j, 0]
                                    if rx > half: rx -= Lbox
                                    elif rx < -half: rx += Lbox
                                    ry = pos[i, 1] - pos[j, 1]
                                    if ry > half: ry -= Lbox
                                    elif ry < -half: ry += Lbox
                                    r2 = rx * rx + ry * ry
                                    if r2 < 1.0:
                                        r = np.sqrt(r2); ssum = svec[i] + svec[j]
                                        if r < ssum and r > 1e-9:
                                            fm = (ssum - r) / r
                                            F[i, 0] += fm * rx; F[i, 1] += fm * ry
                                            F[j, 0] -= fm * rx; F[j, 1] -= fm * ry
                                            # second nonreciprocal channel (unit chi2)
                                            h = (svec[i] - svec[j]) / ssum
                                            f2c = h * fm
                                            F2[i, 0] += f2c * rx; F2[i, 1] += f2c * ry
                                            F2[j, 0] += f2c * rx; F2[j, 1] += f2c * ry
                                            F[i, 0] += chi2 * f2c * rx; F[i, 1] += chi2 * f2c * ry
                                            F[j, 0] += chi2 * f2c * rx; F[j, 1] += chi2 * f2c * ry
                                        t = r2 + l2[j]
                                        gi = alpha * l4[j] / (t * t * np.sqrt(t))
                                        t = r2 + l2[i]
                                        gj = alpha * l4[i] / (t * t * np.sqrt(t))
                                        gbar = 0.5 * (gi + gj)
                                        F[i, 0] -= gbar * rx; F[i, 1] -= gbar * ry
                                        F[j, 0] += gbar * rx; F[j, 1] += gbar * ry
                                        # first nonreciprocal channel (unit chi)
                                        f1c = -(gi - gbar)
                                        F1[i, 0] += f1c * rx; F1[i, 1] += f1c * ry
                                        F1[j, 0] += f1c * rx; F1[j, 1] += f1c * ry
                                        F[i, 0] += chi * f1c * rx; F[i, 1] += chi * f1c * ry
                                        F[j, 0] += chi * f1c * rx; F[j, 1] += chi * f1c * ry
                                j = nxt[j]
                    i = nxt[i]
        for i in range(N):
            inv = 1.0 / svec[i]
            ddx = dt * inv * F[i, 0] + nstd[i] * np.random.randn()
            ddy = dt * inv * F[i, 1] + nstd[i] * np.random.randn()
            accJ1 += F1[i, 0] * ddx + F1[i, 1] * ddy
            accJ2 += F2[i, 0] * ddx + F2[i, 1] * ddy
            pos[i, 0] = (pos[i, 0] + ddx) % Lbox
            pos[i, 1] = (pos[i, 1] + ddy) % Lbox
    T = nsteps * dt
    return accJ1 / T, accJ2 / T


def job(a):
    z = np.load(a["src"])
    # SHORT windows from MANY independent starts. Common random numbers only reduce variance
    # while the paired trajectories stay close; in a chaotic system a long window lets them
    # decorrelate completely and the pairing buys nothing. At T=5e4 the cross-coefficient L21
    # was unresolvable for exactly this reason. Short windows keep the pair correlated, and
    # independent snapshots supply the averaging instead.
    pos = np.ascontiguousarray(z["snaps"][a.get("snap", -1)]).copy()
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    nsteps = int(round(a["T"] / a["dt"]))
    J1, J2 = _run2(pos, svec, lvec**2, lvec**4, nsteps, a["dt"], L,
                   float(z["alpha"]), a["chi"], a["chi2"], SIGMA, a["seed"])
    return dict(chi=a["chi"], chi2=a["chi2"], seed=a["seed"], J1=J1, J2=J2, N=len(svec))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--chi0", type=float, default=0.5)
    p.add_argument("--chi20", type=float, default=0.5)
    p.add_argument("--delta", type=float, default=0.25)
    p.add_argument("--T", type=float, default=5e4)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--pattern", default="data/*_x1_*N1000_*seed[12345].npz")
    p.add_argument("--nstart", type=int, default=0,
                   help="if >0, use this many short windows from independent snapshots")
    p.add_argument("--workers", type=int, default=8)
    a = p.parse_args()
    srcs = sorted(glob.glob(a.pattern))[: (a.nstart or a.seeds)]
    if not srcs:
        raise SystemExit("no equilibrated N=1000 configurations found in data/")
    d = a.delta
    pts = [(a.chi0, a.chi20 - d, "J1-"), (a.chi0, a.chi20 + d, "J1+"),
           (a.chi0 - d, a.chi20, "J2-"), (a.chi0 + d, a.chi20, "J2+")]
    if a.nstart:
        nsnap_each = max(1, a.nstart // max(1, len(srcs)))
        starts = [(s, si) for s in srcs
                  for si in ([0] if nsnap_each == 1 else
                             np.linspace(50, 99, nsnap_each, dtype=int))]
        jobs = [dict(src=s, snap=int(si), chi=c, chi2=c2, seed=2000 + k, T=a.T, dt=a.dt)
                for k, (s, si) in enumerate(starts) for (c, c2, _) in pts]
    else:
        jobs = [dict(src=s, chi=c, chi2=c2, seed=1000 + k, T=a.T, dt=a.dt)
                for k, s in enumerate(srcs) for (c, c2, _) in pts]
    print(f"{len(jobs)} runs ({len(pts)} grid points x {len(srcs)} seeds), {a.workers} workers")
    print(f"operating point (chi, chi2) = ({a.chi0}, {a.chi20}), delta = {d}, T = {a.T:g}")
    with mp.Pool(a.workers) as pool:
        R = pool.map(job, jobs)
    import pandas as pd
    df = pd.DataFrame(R); df.to_csv("onsager.csv", index=False)
    df = df.sort_values("seed")
    key = lambda c, c2: (abs(df.chi - c) < 1e-9) & (abs(df.chi2 - c2) < 1e-9)
    L12 = (df[key(a.chi0, a.chi20 + d)].sort_values("seed").J1.values
           - df[key(a.chi0, a.chi20 - d)].sort_values("seed").J1.values) / (2 * d)
    L21 = (df[key(a.chi0 + d, a.chi20)].sort_values("seed").J2.values
           - df[key(a.chi0 - d, a.chi20)].sort_values("seed").J2.values) / (2 * d)
    print("\n=== Onsager cross-coefficients (paired, common random numbers) ===")
    print(f"  L12 = dJ1/dchi2 = {L12.mean():+.5e} +- {L12.std(ddof=1)/np.sqrt(len(L12)):.2e}")
    print(f"  L21 = dJ2/dchi  = {L21.mean():+.5e} +- {L21.std(ddof=1)/np.sqrt(len(L21)):.2e}")
    dif = L12 - L21
    se = dif.std(ddof=1) / np.sqrt(len(dif))
    print(f"  difference      = {dif.mean():+.5e} +- {se:.2e}   t = {dif.mean()/se:+.2f}")
    rel = abs(dif.mean()) / (0.5 * abs(L12.mean() + L21.mean()))
    print(f"  |L12-L21| / mean(|L|) = {rel:.3f}")
    print(f"\n  {'RECIPROCITY HOLDS (consistent with Onsager, tier II)' if abs(dif.mean()/se) < 2 else 'RECIPROCITY VIOLATED -- beyond linear response'}")
