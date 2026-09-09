"""T_eff with and without whole-system drift, from the SAME starts.

Nonreciprocal forces violate Newton's third law, so the system carries net momentum and its
centre of mass translates ballistically. The response in fdt.py is a paired difference under
common random numbers, so that drift cancels there -- but it does NOT cancel in the fluctuation,
which is a single unperturbed realisation. Pairing a drift-contaminated fluctuation with a
drift-free response inflates T_eff and makes it grow with the window.

This measures both consistently: fluctuation of the species coordinate relative to the
system centre of mass, against the same paired response.
"""
import argparse, glob
import numpy as np
import multiprocessing as mp
from sim.simulate import SIGMA
from analysis.fdt import _run_fdt

kT = SIGMA**2 / 2


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
                      alpha, a["chi"], SIGMA, a["seed"])
        nsp = int(sel.sum())
        out[f"mu_{lbl}"] = (d1[sel, 0] - d0[sel, 0]).sum() / (a["f"] * nsp)
        raw = d0[sel].sum(axis=0)
        rel = (d0[sel] - d0.mean(axis=0)).sum(axis=0)
        out[f"vraw_{lbl}"] = 0.5 * (raw[0]**2 + raw[1]**2) / nsp
        out[f"vrel_{lbl}"] = 0.5 * (rel[0]**2 + rel[1]**2) / nsp
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--chi", type=float, default=1.5)
    p.add_argument("--f", type=float, default=2e-4)
    p.add_argument("--nstart", type=int, default=150)
    p.add_argument("--workers", type=int, default=11)
    p.add_argument("--dt", type=float, default=0.05)
    a = p.parse_args()
    srcs = sorted(glob.glob(f"data/*_x{a.chi:g}_*N1000_*seed[12345].npz"))
    per = max(1, a.nstart // max(1, len(srcs)))
    starts = [(s, int(si)) for s in srcs for si in np.linspace(50, 99, per, dtype=int)]
    print(f"chi={a.chi:g}   {len(starts)} starts   kT={kT:.3e}")
    print("    T     species   T_eff/kT (raw)   T_eff/kT (drift-removed)")
    res = {}
    for T in (200.0, 400.0, 800.0):
        jobs = [dict(src=s, snap=si, chi=a.chi, f=a.f, T=T, dt=a.dt, seed=8000 + k)
                for k, (s, si) in enumerate(starts)]
        with mp.Pool(a.workers) as pool:
            R = pool.map(job, jobs)
        res[T] = {}
        for lbl in ("L", "S"):
            mu = np.array([r[f"mu_{lbl}"] for r in R])
            vr = np.array([r[f"vraw_{lbl}"] for r in R])
            vl = np.array([r[f"vrel_{lbl}"] for r in R])
            traw, trel = vr / (2 * mu), vl / (2 * mu)
            res[T][lbl] = (traw.mean(), trel.mean(), vr.mean(), vl.mean())
            print(f"  {T:6.0f}      {lbl}        {traw.mean()/kT:6.2f} "
                  f"+- {traw.std(ddof=1)/np.sqrt(len(traw))/kT:.2f}"
                  f"        {trel.mean()/kT:6.2f} +- {trel.std(ddof=1)/np.sqrt(len(trel))/kT:.2f}",
                  flush=True)
    lt = np.log([200.0, 400.0, 800.0])
    print("\n  exponents over the window:")
    for lbl in ("L", "S"):
        for k, nm in ((2, "<dX^2> raw"), (3, "<dX^2> drift-removed"), (0, "T_eff raw"), (1, "T_eff drift-removed")):
            al = np.polyfit(lt, np.log([res[T][lbl][k] for T in (200.0, 400.0, 800.0)]), 1)[0]
            print(f"    {lbl}  {nm:24s} alpha = {al:+.3f}")
