"""Is the reciprocal (chi=0) state kinetically arrested, or is it the true equilibrium?

At paper scale chi=0 plateaus at largest-cluster-fraction 0.30 while chi=0.25 condenses to
0.91. Two readings:
  (A) KINETIC ARREST -- the condensed state is also equilibrium-stable, but chi=0 cannot
      reach it because large clusters diffuse too slowly to coalesce. Weak activity unjams it.
  (B) GENUINELY NONEQUILIBRIUM -- the dispersed state is the true chi=0 equilibrium, and
      chi=0.25 condensation is an activity-sustained structure.

These differ in a decisive way: take the CONDENSED configuration and evolve it at chi=0.
Under (A) it stays condensed; under (B) it redisperses toward lcf ~ 0.30.

  python arrest_test.py --workers 8
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp, os
import numpy as np, pandas as pd
import networkx as nx
from sim.simulate import _run, SIGMA
from analysis.analyze import contact_graph


def structure(pos, svec, L):
    G = contact_graph(pos, svec, L)
    comps = [c for c in nx.connected_components(G) if len(c) >= 2]
    return len(comps), max((len(c) for c in comps), default=0) / len(svec)


def job(a):
    src, chi_new, T, tag = a["src"], a["chi"], a["T"], a["tag"]
    z = np.load(src)
    pos = np.ascontiguousarray(z["snaps"][-1]).copy()
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    seed = int(z["seed"])
    n0, l0 = structure(pos, svec, L)
    dt = 0.05
    ns = int(round(T / dt)); se = max(1, ns // 20); nk = ns // se + 1
    snaps = np.zeros((nk, len(svec), 2)); snaps[0] = pos; heat = np.zeros(nk)
    _run(pos, svec, lvec, lvec**2, lvec**4, ns, dt, L, float(z["alpha"]), chi_new,
         SIGMA, seed + 4242, se, snaps, heat)
    traj = [structure(snaps[k], svec, L) for k in range(0, nk, max(1, nk // 5))]
    n1, l1 = structure(pos, svec, L)
    return dict(tag=tag, seed=seed, chi_new=chi_new, n_cl0=n0, lcf0=l0, n_cl1=n1, lcf1=l1,
                lcf_traj=";".join(f"{l:.3f}" for _, l in traj))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--T", type=float, default=6e4)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--datadir", default="data_paper")
    a = p.parse_args()
    disp = sorted(glob.glob(os.path.join(a.datadir, "*_x0_*.npz")))      # dispersed, lcf~0.30
    cond = sorted(glob.glob(os.path.join(a.datadir, "*_x0.25_*.npz")))   # condensed, lcf~0.91
    jobs = []
    for f in cond:
        jobs.append(dict(src=f, chi=0.0,  T=a.T, tag="condensed -> chi=0   (KEY)"))
        jobs.append(dict(src=f, chi=0.25, T=a.T, tag="condensed -> chi=0.25 (ctrl)"))
    for f in disp:
        jobs.append(dict(src=f, chi=0.0,  T=a.T, tag="dispersed -> chi=0    (ctrl)"))
        jobs.append(dict(src=f, chi=0.25, T=a.T, tag="dispersed -> chi=0.25 (KEY)"))
    print(f"{len(jobs)} jobs, T={a.T:g}")
    with mp.Pool(a.workers) as pool:
        rows = pool.map(job, jobs)
    df = pd.DataFrame(rows); df.to_csv("arrest_test.csv", index=False)
    print("\n  transition                      lcf: start -> end     n_cl: start -> end")
    for tag, d in df.groupby("tag"):
        print(f"  {tag:32s} {d.lcf0.mean():.3f} -> {d.lcf1.mean():.3f}      "
              f"{d.n_cl0.mean():6.1f} -> {d.n_cl1.mean():6.1f}")
    print("\n  lcf trajectories (condensed -> chi=0):")
    for _, r in df[df.tag.str.startswith("condensed -> chi=0 ")].iterrows():
        print(f"    seed{r.seed}: {r.lcf_traj}")


if __name__ == "__main__":
    main()
