"""Aspect-ratio control: is the square-box approximation in the source replicate justified?

The source domain is 648 x 360 um (aspect ratio 1.8); our source-specification replicate used an
equal-area square. This runs the same N, density and duration as the converged square-box case at
N=4000, in a 1.8:1 rectangle of identical area, so the comparison isolates geometry alone.
"""
from __future__ import annotations
import math, multiprocessing as mp, os
import numpy as np
from sim.simulate import setup, SIGMA
from simulate_rect import run_rect

AREA, AR = 24.0 ** 2, 648 / 360
LY = math.sqrt(AREA / AR)
LX = AR * LY


def job(a):
    tag = f"rect_x{a['chi']:g}_N{a['N']}_AR{AR:.1f}_T{a['T']:g}_seed{a['seed']}"
    out = os.path.join("data_rect", tag + ".npz")
    if os.path.exists(out):
        return f"skip {tag}"
    rng = np.random.default_rng(a["seed"])
    types, svec, lvec = setup("bidisperse", rng, a["N"], 1 / 1.5)
    pos = np.column_stack([rng.uniform(0, LX, a["N"]), rng.uniform(0, LY, a["N"])])
    ns = int(round(a["T"] / a["dt"]))
    se = max(1, ns // a["nsnap"])
    snaps = np.zeros((ns // se + 1, a["N"], 2)); snaps[0] = pos
    run_rect(pos, svec, lvec, ns, a["dt"], LX, LY, 0.005, a["chi"], SIGMA,
             a["seed"], se, snaps)
    np.savez_compressed(out, snaps=snaps, svec=svec, lvec=lvec, types=types,
                        dt_snap=se * a["dt"], Lbox=LX, Lx=LX, Ly=LY, alpha=0.005,
                        s_ratio=1 / 1.5, chi=a["chi"], dt=a["dt"], seed=a["seed"],
                        case="bidisperse")
    return f"done {tag}"


if __name__ == "__main__":
    os.makedirs("data_rect", exist_ok=True)
    print(f"rectangle {LX:.2f} x {LY:.2f} (area {LX*LY:.0f}, AR {AR:.2f}) "
          f"vs converged square 24.00 x 24.00 (area {AREA:.0f})")
    jobs = [dict(N=4000, chi=1.5, T=8e5, dt=0.05, nsnap=200, seed=s) for s in (1, 2, 3)]
    with mp.Pool(3) as pool:
        for m in pool.imap_unordered(job, jobs):
            print(" ", m, flush=True)
