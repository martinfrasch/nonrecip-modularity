"""Committor of fragment clusters (OPEN_QUESTIONS_PLAN.md, track C3).

From stored dense-snapshot configurations, relaunch the full system with fresh noise M times for
a window T and follow every fragment (size ≥ 3, not the condensate) present at t = 0 by
plurality tracking. Its fate is the first of
  grow      tracked size ≥ S_big
  shrink    tracked size ≤ S_small
  absorbed  its plurality successor is the condensate
  open      none within T
q(S) = P(grow ∪ absorbed) / P(grow ∪ absorbed ∪ shrink) is the committor to the condensed side;
q_free(S) excludes absorption. A critical nucleus needs q(S*) ≈ ½ with q monotone in S.

  python committor.py --pattern "data_paper/*cont2_dense500.npz" --starts 4 --M 25 --T 500 --workers 4
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp, os
import numpy as np
import pandas as pd
from simulate import _run, SIGMA
from cluster_kinetics import labels

S_SMALL, S_BIG = 2, 20


def one(a):
    f, t0, m, T, dt_snap = a["file"], a["t0"], a["m"], a["T"], a["dt_snap"]
    z = np.load(f)
    sv, lv, L, al = z["svec"], z["lvec"], float(z["Lbox"]), float(z["alpha"])
    chi = float(z["chi"]); dt = float(z["dt"]); seed = int(z["seed"])
    pos = np.ascontiguousarray(z["snaps"][t0]).copy()
    ns = int(T / dt); se = int(dt_snap / dt)
    snaps = np.zeros((ns // se + 1, len(sv), 2)); heat = np.zeros(len(snaps))
    snaps[0] = pos
    _run(pos, sv, lv, lv**2, lv**4, ns, dt, L, al, chi, SIGMA, 100000 + 1000 * seed + 37 * t0 + m, se, snaps, heat)
    lab0, sz0 = labels(snaps[0], sv, L)
    frag = [k for k in range(1, len(sz0)) if sz0[k] >= 3]
    members = {k: np.where(lab0 == k)[0] for k in frag}
    fate = {k: None for k in frag}; when = {k: None for k in frag}; size_t = {k: sz0[k] for k in frag}
    for s in range(1, len(snaps)):
        lab, sz = labels(snaps[s], sv, L)
        for k in frag:
            if fate[k] is not None:
                continue
            succ, cnt = np.unique(lab[members[k]], return_counts=True)
            main = succ[np.argmax(cnt)]
            S = sz[main]
            if main == 0:
                fate[k], when[k] = "absorbed", s
            elif S >= S_BIG:
                fate[k], when[k] = "grow", s
            elif S <= S_SMALL:
                fate[k], when[k] = "shrink", s
            else:
                members[k] = np.where(lab == main)[0]   # follow the plurality successor
    return [dict(file=os.path.basename(f), chi=chi, seed=seed, t0=t0, m=m, S0=int(sz0[k]),
                 fate=fate[k] or "open", when=(when[k] or 0) * dt_snap) for k in frag]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--starts", type=int, default=4)
    p.add_argument("--M", type=int, default=25)
    p.add_argument("--T", type=float, default=500.0)
    p.add_argument("--dt-snap", type=float, default=50.0)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--out", default="committor.csv")
    a = p.parse_args()
    jobs = []
    for f in sorted(set(x for pat in a.pattern for x in glob.glob(pat))):
        n = len(np.load(f)["snaps"])
        for t0 in np.linspace(n // 2, n - 1, a.starts, dtype=int):
            for m in range(a.M):
                jobs.append(dict(file=f, t0=int(t0), m=m, T=a.T, dt_snap=a.dt_snap))
    print(f"{len(jobs)} launches", flush=True)
    with mp.Pool(a.workers) as pool:
        for k, rows in enumerate(pool.imap_unordered(one, jobs)):
            pd.DataFrame(rows).to_csv(a.out, mode="a", header=not os.path.exists(a.out), index=False)
            if k % 10 == 0:
                print(f"{k + 1}/{len(jobs)} done", flush=True)


if __name__ == "__main__":
    main()
