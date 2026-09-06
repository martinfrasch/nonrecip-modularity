"""
Agent-based model of Hara et al., PRL 137, 068302 (2026) / arXiv:2509.23164,
End Matter Sec. II, Eqs. (7)-(9). Nondimensional units: lambda = 9 um, tau = 0.01 s.

  dx_i/dt = xi_i(t) + (1/s_i) [ sum_j F^col_ij + sum_j F^EHD_ij ]

  F^col_ij : reciprocal soft-core linear spring, natural length s_i + s_j
  F^EHD_ij = -alpha * l_j^4 / (r^2 + l_j^2)^{5/2} * r_vec   (cutoff r <= 1)
             NONRECIPROCAL: force on i from j scales with l_j (the *other*
             particle's EHD radius), so unequal particles break action-reaction
             symmetry. The antisymmetric part of this coupling is the object
             of interest for the NWAP directed-graph extension.

  noise: <xi_i xi_j> = (sigma^2 / s_i) delta_ij delta(t-t'), sigma = 2.6e-3

Paper baseline: alpha = 0.005, l_L = s_L = 1/6 (1.5 um), l_S = s_S = 1/9 (1 um),
steric polydispersity eps ~ N(0, 1/30), head-large mix N_L : N_S = 1 : 3.

Usage:
  python simulate.py --case bidisperse --seed 12
  python simulate.py --case monodisperse --seed 11
  python simulate.py --sweep baseline          # 5 seeds x both cases
  python simulate.py --sweep alpha             # alpha in {0.003,0.005,0.010,0.015}
  python simulate.py --sweep sizeratio         # s_II/s_I sweep at fixed EHD radii
                                               # (head-large -> tail-large, cf. Fig 4 / S5)
  python simulate.py --case bidisperse --N 4000 --box 24 --T 4e5   # paper-scale run
"""
from __future__ import annotations
import argparse
import itertools
import multiprocessing as mp
import os
import time

import numpy as np
from numba import njit

RCUT = 1.0
SIGMA = 2.6e-3
DEFAULT_DT = 0.05          # conservative on a workstation; container run used 0.12


@njit(cache=True, fastmath=True)
def _run(pos, svec, lvec, l2, l4, nsteps, dt, Lbox, alpha, chi, sigma, seed, snap_every, snaps, heat):
    np.random.seed(seed)
    N = pos.shape[0]
    ncell = int(Lbox / RCUT)
    csize = Lbox / ncell
    head = np.full(ncell * ncell, -1, np.int64)
    nxt = np.full(N, -1, np.int64)
    F = np.zeros((N, 2))
    half = 0.5 * Lbox
    nstd = np.empty(N)
    for i in range(N):
        nstd[i] = sigma * np.sqrt(dt / svec[i])
    isnap = 1
    # Stratonovich heat  sum_i F_i o dx_i  -- the entropy production / current probe.
    # Exactly zero (in the mean) under detailed balance; positive iff the dynamics
    # break time-reversal. Midpoint rule needs F at both ends of a step, so the
    # accumulation runs one step behind.
    Fprev = np.zeros((N, 2))
    dxs = np.zeros((N, 2))
    acc = 0.0
    for step in range(1, nsteps + 1):
        head[:] = -1
        for i in range(N):
            c = (int(pos[i, 0] / csize) % ncell) * ncell + (int(pos[i, 1] / csize) % ncell)
            nxt[i] = head[c]
            head[c] = i
        F[:] = 0.0
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
                                    if rx > half:
                                        rx -= Lbox
                                    elif rx < -half:
                                        rx += Lbox
                                    ry = pos[i, 1] - pos[j, 1]
                                    if ry > half:
                                        ry -= Lbox
                                    elif ry < -half:
                                        ry += Lbox
                                    r2 = rx * rx + ry * ry
                                    if r2 < 1.0:
                                        r = np.sqrt(r2)
                                        ssum = svec[i] + svec[j]
                                        if r < ssum and r > 1e-9:
                                            fm = (ssum - r) / r
                                            F[i, 0] += fm * rx
                                            F[i, 1] += fm * ry
                                            F[j, 0] -= fm * rx
                                            F[j, 1] -= fm * ry
                                        # INDEX CONVENTION: gi is the coefficient acting ON i,
                                        # and it is built from l[j] -- the OTHER particle's EHD
                                        # radius. That asymmetry IS the nonreciprocity of the
                                        # source model; do not "fix" it to l[i].
                                        t = r2 + l2[j]
                                        gi = alpha * l4[j] / (t * t * np.sqrt(t))   # acts on i
                                        t = r2 + l2[i]
                                        gj = alpha * l4[i] / (t * t * np.sqrt(t))   # acts on j
                                        # reciprocity mixing: chi=1 -> original (Hara et al.),
                                        # chi=0 -> exactly reciprocal. The symmetric part
                                        # (gi+gj)/2 is invariant in chi; only the
                                        # antisymmetric part scales, linearly.
                                        gbar = 0.5 * (gi + gj)
                                        ci = gbar + chi * (gi - gbar)
                                        cj = gbar + chi * (gj - gbar)
                                        F[i, 0] -= ci * rx
                                        F[i, 1] -= ci * ry
                                        F[j, 0] += cj * rx
                                        F[j, 1] += cj * ry
                                j = nxt[j]
                    i = nxt[i]
        if step > 1:
            for i in range(N):
                acc += 0.5 * ((Fprev[i, 0] + F[i, 0]) * dxs[i, 0]
                              + (Fprev[i, 1] + F[i, 1]) * dxs[i, 1])
        for i in range(N):
            Fprev[i, 0] = F[i, 0]
            Fprev[i, 1] = F[i, 1]
        for i in range(N):
            inv = 1.0 / svec[i]
            ddx = dt * inv * F[i, 0] + nstd[i] * np.random.randn()
            ddy = dt * inv * F[i, 1] + nstd[i] * np.random.randn()
            dxs[i, 0] = ddx
            dxs[i, 1] = ddy
            pos[i, 0] = (pos[i, 0] + ddx) % Lbox
            pos[i, 1] = (pos[i, 1] + ddy) % Lbox
        if step % snap_every == 0:
            snaps[isnap] = pos
            heat[isnap] = acc
            isnap += 1
    return pos


def setup(case: str, rng: np.random.Generator, N: int, s_ratio_II_over_I: float = 1 / 1.5,
          frac_large: float = 0.25):
    """Build particle arrays.

    case 'bidisperse': type I (large, l=s=1/6) : type II (small) = 1 : 3, head-large
                       by default (s_II/s_I = 0.667). Set s_ratio_II_over_I = 1.5
                       for the tail-large control of Fig. 4(b) (steric radii swapped,
                       EHD radii fixed).
    case 'monodisperse': all L particles at matched area packing (~30%).
    """
    if case == "bidisperse":
        NL = int(round(N * frac_large))     # paper: 5000/22000 = 0.227; our earlier runs used 0.25
        NS = N - NL
        types = np.array([0] * NL + [1] * NS)
        sI = 1 / 6
        sII = sI * s_ratio_II_over_I
        s0 = np.where(types == 0, sI, sII)
        l0 = np.where(types == 0, 1 / 6, 1 / 9)     # EHD radii fixed (paper S5/S8 protocol)
    elif case == "monodisperse":
        types = np.zeros(N, np.int64)
        s0 = np.full(N, 1 / 6)
        l0 = np.full(N, 1 / 6)
    else:
        raise ValueError(case)
    svec = s0 * (1 + rng.normal(0, 1 / 30, len(types)))
    lvec = l0.astype(float)
    return types, svec, lvec


def run_one(job: dict) -> str:
    case, seed = job["case"], job["seed"]
    N, Lbox, T, dt, alpha, nsnap = job["N"], job["box"], job["T"], job["dt"], job["alpha"], job["nsnap"]
    sr = job.get("s_ratio", 1 / 1.5)
    chi = job.get("chi", 1.0)
    fl = job.get("frac_large", 0.25)
    outdir = job["outdir"]
    # frac_large appears in the tag only when non-default, so earlier runs stay resumable
    fltag = "" if abs(fl - 0.25) < 1e-9 else f"_fl{fl:g}"
    tag = f"{case}_a{alpha:g}_sr{sr:.3f}_x{chi:g}{fltag}_N{N}_L{Lbox:g}_T{T:g}_seed{seed}"
    out = os.path.join(outdir, tag + ".npz")
    if os.path.exists(out):
        return f"skip {tag} (exists)"
    rng = np.random.default_rng(seed)
    # monodisperse at matched packing: half the particle count of the bidisperse mix
    Nc = N if case == "bidisperse" else N // 2
    types, svec, lvec = setup(case, rng, Nc, sr, fl)
    pos = rng.uniform(0, Lbox, (Nc, 2))
    nsteps = int(round(T / dt))
    snap_every = max(1, nsteps // nsnap)
    nkeep = nsteps // snap_every + 1
    snaps = np.zeros((nkeep, Nc, 2))
    snaps[0] = pos
    t0 = time.time()
    heat = np.zeros(nkeep)
    _run(pos, svec, lvec, lvec**2, lvec**4, nsteps, dt, Lbox, alpha, chi, SIGMA, seed,
         snap_every, snaps, heat)
    np.savez_compressed(out, snaps=snaps, heat=heat, svec=svec, lvec=lvec, types=types,
                        dt_snap=snap_every * dt, Lbox=Lbox, alpha=alpha, s_ratio=sr,
                        chi=chi, frac_large=fl, dt=dt, seed=seed, case=case)
    return f"done {tag}  wall={time.time() - t0:.0f}s"


def build_jobs(args) -> list[dict]:
    base = dict(N=args.N, box=args.box, T=args.T, dt=args.dt, alpha=args.alpha,
                chi=args.chi, frac_large=args.frac_large, nsnap=args.nsnap, outdir=args.outdir)
    if args.sweep == "baseline":
        return [dict(base, case=c, seed=s)
                for c, s in itertools.product(["monodisperse", "bidisperse"], range(1, 6))]
    if args.sweep == "alpha":
        return [dict(base, case="bidisperse", seed=s, alpha=a)
                for a, s in itertools.product([0.003, 0.005, 0.010, 0.015], range(1, 4))]
    if args.sweep == "chi":
        # THE reciprocity sweep: chi scales the antisymmetric coupling alone,
        # holding the symmetric part identical. chi=0 is a same-particle
        # reciprocal control; chi=1 reproduces baseline bidisperse.
        return [dict(base, case="bidisperse", seed=s, chi=x)
                for x, s in itertools.product([0.0, 0.25, 0.5, 0.75, 1.0, 1.5], range(1, 6))]
    if args.sweep == "validate":
        # replicate at the paper's own specification (End Matter / SI Sec. III)
        return [dict(base, case="bidisperse", seed=s, chi=1.0) for s in range(1, args.seeds + 1)]
    if args.sweep == "chibox":
        # box-scaling at fixed density: N/L^2 must match across boxes.
        # If a finite characteristic cluster size S* exists, largest-cluster fraction
        # falls as S*/N ~ 1/L^2. If the system truly phase-separates, lcf is L-independent.
        chis = [float(x) for x in args.chi_list.split(",")]
        return [dict(base, case="bidisperse", seed=s, chi=c)
                for c, s in itertools.product(chis, range(1, args.seeds + 1))]
    if args.sweep == "chipaper":
        # paper-scale reciprocity test: full chi curve at 3 seeds, plus a
        # monodisperse reference at matched scale (mono runs N/2 particles)
        jobs = [dict(base, case="bidisperse", seed=s, chi=x)
                for x, s in itertools.product([0.0, 0.25, 0.5, 0.75, 1.0, 1.5], range(1, 4))]
        jobs += [dict(base, case="monodisperse", seed=s, chi=1.0) for s in range(1, 4)]
        return jobs
    if args.sweep == "sizeratio":
        # 0.667 = head-large (fragmenting), 1.0 = symmetric steric, 1.5 = tail-large (aggregating)
        return [dict(base, case="bidisperse", seed=s, s_ratio=r)
                for r, s in itertools.product([1 / 1.5, 1 / 1.25, 1.0, 1.25, 1.5], range(1, 4))]
    return [dict(base, case=args.case, seed=args.seed)]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", choices=["monodisperse", "bidisperse"], default="bidisperse")
    p.add_argument("--seed", type=int, default=12)
    p.add_argument("--sweep", choices=["baseline", "alpha", "sizeratio", "chi", "chipaper", "chibox", "validate"], default=None)
    p.add_argument("--N", type=int, default=1000, help="bidisperse particle count (mono uses N/2); paper scale ~4000")
    p.add_argument("--box", type=float, default=12.0, help="domain edge in units of lambda=9um; paper Fig 4 uses 24")
    p.add_argument("--T", type=float, default=1e5, help="total nondim time (1e5 = 1000 s); paper steady-state stats need >=4e5")
    p.add_argument("--dt", type=float, default=DEFAULT_DT)
    p.add_argument("--alpha", type=float, default=0.005)
    p.add_argument("--chi", type=float, default=1.0,
                   help="reciprocity mixing: 0 = exactly reciprocal, 1 = Hara et al. force")
    p.add_argument("--nsnap", type=int, default=100)
    p.add_argument("--frac-large", type=float, default=0.25,
                   help="fraction of type-I (large) particles; paper uses 5000/22000 = 0.2273")
    p.add_argument("--chi-list", default="1.0,1.5", help="chi values for --sweep chibox")
    p.add_argument("--seeds", type=int, default=2, help="seeds per level for --sweep chibox")
    p.add_argument("--outdir", default="data")
    p.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 1))
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    jobs = build_jobs(args)
    print(f"{len(jobs)} job(s), {args.workers} worker(s)")
    if len(jobs) == 1 or args.workers == 1:
        for j in jobs:
            print(run_one(j), flush=True)
    else:
        with mp.Pool(args.workers) as pool:
            for msg in pool.imap_unordered(run_one, jobs):
                print(msg, flush=True)
