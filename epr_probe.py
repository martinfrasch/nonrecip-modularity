"""Entropy production rate: the direct test of whether the dynamics break time-reversal.

EPR is exactly zero under detailed balance and positive iff a probability current flows,
so it is the measurement that H_B ("the signature appears in probability currents") needs
and that edge turnover / sigma^2_v cannot provide -- both are nonzero in the provably
current-free chi=0 system.

The Stratonovich heat estimator subtracts two nearly-cancelling O(mu|F|^2) terms, so its
discretisation residual swamps the signal at the production dt=0.05 (see AUDIT.md). It is
usable only at much smaller dt. EPR is an intensive steady-state quantity, so rather than
re-running long trajectories we restart from the already-equilibrated final snapshot of an
existing production run and integrate a short window at small dt.

  python epr_probe.py --dt 6.25e-4 --T 600 --datadir data
"""
from __future__ import annotations
import argparse, glob, os, re
import numpy as np
import pandas as pd
from simulate import _run, SIGMA

T_EFF = SIGMA**2 / 2          # D_i/mu_i = sigma^2/2, uniform -> effective temperature


def probe(path, dt, T, nwin=20):
    d = np.load(path)
    pos = np.ascontiguousarray(d["snaps"][-1]).copy()   # equilibrated config
    svec, lvec, L = d["svec"], d["lvec"], float(d["Lbox"])
    alpha = float(d["alpha"])
    chi = float(d["chi"]) if "chi" in d.files else 1.0
    seed = int(d["seed"])
    N = len(svec)
    nsteps = int(round(T / dt))
    snap_every = max(1, nsteps // nwin)
    nkeep = nsteps // snap_every + 1
    snaps = np.zeros((nkeep, N, 2)); snaps[0] = pos
    heat = np.zeros(nkeep)
    _run(pos, svec, lvec, lvec**2, lvec**4, nsteps, dt, L, alpha, chi,
         SIGMA, seed + 9000, snap_every, snaps, heat)
    lo = nkeep // 2                                     # drop the dt-switch transient
    dQ = heat[-1] - heat[lo]
    dtime = (nkeep - 1 - lo) * snap_every * dt
    return dict(file=os.path.basename(path), chi=chi, seed=seed, N=N,
                epr=dQ / dtime / T_EFF / N)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--datadir", default="data")
    p.add_argument("--dt", type=float, default=6.25e-4)
    p.add_argument("--T", type=float, default=600.0)
    p.add_argument("--pattern", default="*_x*_N1000_*.npz")
    p.add_argument("--out", default="epr.csv")
    a = p.parse_args()
    files = sorted(glob.glob(os.path.join(a.datadir, a.pattern)))
    if not files:
        raise SystemExit(f"no files matching {a.pattern} in {a.datadir}")
    rows = []
    for f in files:
        r = probe(f, a.dt, a.T)
        rows.append(r)
        print(f"  chi={r['chi']:<5g} seed={r['seed']}  EPR/particle = {r['epr']:+.4e}", flush=True)
    df = pd.DataFrame(rows); df.to_csv(a.out, index=False)
    g = df.groupby("chi").agg(epr=("epr", "mean"), err=("epr", "sem"), n=("seed", "count"))
    z = g.epr.loc[0.0] if 0.0 in g.index else np.nan
    print(f"\n== EPR per particle (dt={a.dt:g}, T={a.T:g}) ==")
    for chi, r in g.iterrows():
        excess = r.epr - z
        q = f"{excess/chi**2:+.3e}" if chi > 0 else "      --"
        print(f"  chi={chi:<5g} {r.epr:+.4e} +- {r['err']:.2e}   excess over chi=0: {excess:+.3e}   /chi^2: {q}")
    print(f"\n  chi=0 is the exact-zero reference; its value {z:+.3e} is the residual")
    print(f"  discretisation bias. Excess over it is the physical entropy production.")


if __name__ == "__main__":
    main()
