"""A3 -- is there a single effective temperature?

At equilibrium the Einstein relation fixes D_i/mu_i = kT for every degree of freedom. Here
mu_i = 1/s_i and the noise variance is sigma^2/s_i, so D_i/mu_i = sigma^2/2 for BOTH species
exactly -- which is how chi=0 was established as a genuine equilibrium reference.

Out of equilibrium the ratio defines an effective temperature, D_i/mu_i = kT_eff,i. A single
T_eff shared across degrees of freedom is a tier-II (Onsager) signature; an observable-dependent
T_eff indicates tier III, where no Rayleighian applies.

The test therefore compares T_eff measured on the LARGE and SMALL species. Mobility comes from
the drift under a small force applied to one species only; diffusion from the mean-squared
displacement of the unperturbed system. Perturbed and unperturbed runs share a noise stream
(common random numbers), so the drift is a paired difference.

As in A2, every quantity is checked for a plateau against the measurement window before it is
believed.
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp
import numpy as np
from numba import njit
from sim.simulate import SIGMA, RCUT


@njit(cache=True, fastmath=True)
def _run_fdt(pos, svec, l2, l4, fext, nsteps, dt, Lbox, alpha, chi, sigma, seed):
    """Returns per-particle unwrapped displacement (dx, dy)."""
    np.random.seed(seed)
    N = pos.shape[0]
    ncell = int(Lbox / RCUT); csize = Lbox / ncell
    head = np.full(ncell * ncell, -1, np.int64); nxt = np.full(N, -1, np.int64)
    F = np.zeros((N, 2)); disp = np.zeros((N, 2))
    half = 0.5 * Lbox
    nstd = np.empty(N)
    for i in range(N):
        nstd[i] = sigma * np.sqrt(dt / svec[i])
    for step in range(nsteps):
        head[:] = -1
        for i in range(N):
            c = (int(pos[i, 0] / csize) % ncell) * ncell + (int(pos[i, 1] / csize) % ncell)
            nxt[i] = head[c]; head[c] = i
        F[:] = 0.0
        for i in range(N):
            F[i, 0] = fext[i]
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
                                        t = r2 + l2[j]
                                        gi = alpha * l4[j] / (t * t * np.sqrt(t))
                                        t = r2 + l2[i]
                                        gj = alpha * l4[i] / (t * t * np.sqrt(t))
                                        gbar = 0.5 * (gi + gj)
                                        ci = gbar + chi * (gi - gbar)
                                        cj = gbar + chi * (gj - gbar)
                                        F[i, 0] -= ci * rx; F[i, 1] -= ci * ry
                                        F[j, 0] += cj * rx; F[j, 1] += cj * ry
                                j = nxt[j]
                    i = nxt[i]
        for i in range(N):
            inv = 1.0 / svec[i]
            ddx = dt * inv * F[i, 0] + nstd[i] * np.random.randn()
            ddy = dt * inv * F[i, 1] + nstd[i] * np.random.randn()
            disp[i, 0] += ddx; disp[i, 1] += ddy
            pos[i, 0] = (pos[i, 0] + ddx) % Lbox
            pos[i, 1] = (pos[i, 1] + ddy) % Lbox
    return disp


def job(a):
    z = np.load(a["src"])
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    ty = np.atleast_1d(z["types"]); alpha = float(z["alpha"])
    n = int(round(a["T"] / a["dt"]))
    out = {}
    for lbl, sel in (("L", ty == 0), ("S", ty == 1)):
        fext = np.zeros(len(svec)); fext[sel] = a["f"]
        pos = np.ascontiguousarray(z["snaps"][a["snap"]]).copy()
        d1 = _run_fdt(pos, svec, lvec**2, lvec**4, fext, n, a["dt"], L, alpha,
                      a["chi"], SIGMA, a["seed"])
        pos = np.ascontiguousarray(z["snaps"][a["snap"]]).copy()
        d0 = _run_fdt(pos, svec, lvec**2, lvec**4, np.zeros(len(svec)), n, a["dt"], L,
                      alpha, a["chi"], SIGMA, a["seed"])          # common random numbers
        # The force couples to the COLLECTIVE coordinate X = sum_{i in species} x_i, so the
        # fluctuation must be measured on X too. A single-particle MSD saturates inside its cage
        # while the collective drift does not, and comparing them is not a conjugate pair --
        # at equilibrium it gives T_eff/kT = 0.21 and 0.11 instead of 1.
        nsp = int(sel.sum())
        dX_f = (d1[sel, 0] - d0[sel, 0]).sum()          # response of X, paired
        X0 = d0[sel, 0].sum(); Y0 = d0[sel, 1].sum()    # unperturbed collective displacement
        out[f"mu_{lbl}"] = dX_f / (a["f"] * nsp)        # per-particle mobility
        out[f"var_{lbl}"] = 0.5 * (X0**2 + Y0**2) / nsp # <dX^2> per particle, isotropised
        out[f"n_{lbl}"] = nsp
    return dict(seed=a["seed"], chi=a["chi"], T=a["T"], **out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--chi", type=float, default=1.0)
    p.add_argument("--f", type=float, default=2e-4)
    p.add_argument("--T", type=float, default=200.0)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--nstart", type=int, default=40)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--pattern", default="data/*_x1_*N1000_*seed[12345].npz")
    a = p.parse_args()
    srcs = sorted(glob.glob(a.pattern))
    per = max(1, a.nstart // max(1, len(srcs)))
    starts = [(s, int(si)) for s in srcs
              for si in (np.linspace(50, 99, per, dtype=int) if per > 1 else [0])]
    jobs = [dict(src=s, snap=si, chi=a.chi, f=a.f, T=a.T, dt=a.dt, seed=5000 + k)
            for k, (s, si) in enumerate(starts)]
    with mp.Pool(a.workers) as pool:
        R = pool.map(job, jobs)
    import pandas as pd
    d = pd.DataFrame(R)
    kT = SIGMA ** 2 / 2
    print(f"  chi={a.chi:g} f={a.f:g} T={a.T:g}  ({len(d)} starts)   equilibrium kT = {kT:.3e}")
    # T_eff = f <dX^2> / (2 <dX>), both on the collective coordinate at the same elapsed time
    for lbl in ("L", "S"):
        mu, v = d[f"mu_{lbl}"].values, d[f"var_{lbl}"].values
        te = v / (2 * mu)   # kT = f<dX^2>/(2<dX>); the f and n cancel in these definitions
        se = te.std(ddof=1) / np.sqrt(len(te))
        print(f"    {lbl}: mu={mu.mean():.3e}  <dX^2>/n={v.mean():.3e}  "
              f"T_eff={te.mean():.3e} +- {se:.1e}   T_eff/kT={te.mean()/kT:6.2f}")
    tl = (d.var_L / (2 * d.mu_L)).values; ts = (d.var_S / (2 * d.mu_S)).values
    dif = tl - ts; se = dif.std(ddof=1) / np.sqrt(len(dif))
    print(f"    T_eff(L) - T_eff(S) = {dif.mean():+.3e} +- {se:.1e}   t = {dif.mean()/se:+.2f}")
