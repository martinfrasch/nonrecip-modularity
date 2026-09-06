"""Is the collective-coordinate superdiffusion just the condensate translating?

Section 3.10 rests its only converged negative result -- no window-independent effective
temperature -- on <dX^2> ~ t^1.81 for X = sum_{i in species} x_i. But nonreciprocal forces
violate Newton's third law, so the system carries net momentum and the whole aggregate can
translate ballistically. That alone gives <dX^2> ~ t^2 by construction.

Control: recompute the same quantity with the system-wide centre-of-mass displacement removed,
    X_rel = sum_{i in species} (dx_i - <dx>_all)
If the exponent drops toward 1, the superdiffusion is whole-system drift, not anomalous
diffusion of the species coordinate.
"""
import argparse, glob
import numpy as np
import multiprocessing as mp
from simulate import SIGMA
from fdt import _run_fdt


def job(a):
    z = np.load(a["src"])
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    ty = np.atleast_1d(z["types"]); alpha = float(z["alpha"])
    n = int(round(a["T"] / a["dt"]))
    pos = np.ascontiguousarray(z["snaps"][a["snap"]]).copy()
    d = _run_fdt(pos, svec, lvec**2, lvec**4, np.zeros(len(svec)), n, a["dt"], L,
                 alpha, a["chi"], SIGMA, a["seed"])
    out = {}
    dall = d.mean(axis=0)                     # system-wide COM displacement
    for lbl, sel in (("L", ty == 0), ("S", ty == 1)):
        nsp = int(sel.sum())
        raw = d[sel].sum(axis=0)              # species collective displacement
        rel = (d[sel] - dall).sum(axis=0)     # ... with system drift removed
        out[f"raw_{lbl}"] = 0.5 * (raw[0]**2 + raw[1]**2) / nsp
        out[f"rel_{lbl}"] = 0.5 * (rel[0]**2 + rel[1]**2) / nsp
    out["cm2"] = 0.5 * (dall[0]**2 + dall[1]**2) * len(svec)   # system COM, same normalisation
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--chi", type=float, default=1.5)
    p.add_argument("--nstart", type=int, default=100)
    p.add_argument("--workers", type=int, default=11)
    p.add_argument("--dt", type=float, default=0.05)
    a = p.parse_args()
    srcs = sorted(glob.glob(f"data/*_x{a.chi:g}_*N1000_*seed[12345].npz"))
    per = max(1, a.nstart // max(1, len(srcs)))
    starts = [(s, int(si)) for s in srcs for si in np.linspace(50, 99, per, dtype=int)]
    print(f"chi={a.chi:g}  {len(starts)} starts from {len(srcs)} runs")
    Ts = [200.0, 400.0, 800.0]
    res = {}
    for T in Ts:
        jobs = [dict(src=s, snap=si, chi=a.chi, T=T, dt=a.dt, seed=7000 + k)
                for k, (s, si) in enumerate(starts)]
        with mp.Pool(a.workers) as pool:
            R = pool.map(job, jobs)
        res[T] = {k: np.mean([r[k] for r in R]) for k in R[0]}
        r = res[T]
        print(f"  T={T:6.0f}  raw_L={r['raw_L']:.4e}  rel_L={r['rel_L']:.4e}  "
              f"raw_S={r['raw_S']:.4e}  rel_S={r['rel_S']:.4e}  sysCM={r['cm2']:.4e}", flush=True)
    print("\n  exponents  <.> ~ t^alpha  (log-log fit over T = 200, 400, 800)")
    lt = np.log([float(T) for T in Ts])
    for k in ("raw_L", "rel_L", "raw_S", "rel_S", "cm2"):
        al = np.polyfit(lt, np.log([res[T][k] for T in Ts]), 1)[0]
        print(f"    {k:7s} alpha = {al:+.3f}")
    print("\n  raw = species collective coordinate (as used in Sec 3.10)")
    print("  rel = same, with system-wide COM displacement removed")
    print("  cm2 = system COM displacement alone, normalised the same way")
