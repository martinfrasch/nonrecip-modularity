"""Continue an existing run from its final snapshot, to settle whether it had converged.

Writes a new .npz with the same metadata and a _cont{n} tag, so analyze.py treats it as an
independent condition and the two windows can be compared.

  python continue_run.py --pattern "data_paper/*_x0.5_*.npz" --T 4e5 --workers 3
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp, os
import numpy as np
from simulate import _run, SIGMA


def job(a):
    src, T, outdir = a["src"], a["T"], a["outdir"]
    z = np.load(src)
    base = os.path.basename(src)[:-4]
    prev = int(base.split("_cont")[1]) if "_cont" in base else 0
    stem = base.split("_cont")[0]
    out = os.path.join(outdir, f"{stem}_cont{prev+1}.npz")
    if os.path.exists(out):
        return f"skip {os.path.basename(out)}"
    pos = np.ascontiguousarray(z["snaps"][-1]).copy()
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    dt = float(z["dt"]); alpha = float(z["alpha"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    seed = int(z["seed"])
    nsnap = len(z["snaps"]) - 1
    ns = int(round(T / dt)); se = max(1, ns // nsnap); nk = ns // se + 1
    snaps = np.zeros((nk, len(svec), 2)); snaps[0] = pos; heat = np.zeros(nk)
    _run(pos, svec, lvec, lvec**2, lvec**4, ns, dt, L, alpha, chi, SIGMA,
         seed + 7777, se, snaps, heat)
    np.savez_compressed(out, snaps=snaps, heat=heat, svec=svec, lvec=lvec,
                        types=z["types"], dt_snap=se*dt, Lbox=L, alpha=alpha,
                        s_ratio=float(z["s_ratio"]), chi=chi, dt=dt, seed=seed,
                        case=str(z["case"]))
    return f"done {os.path.basename(out)}"


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", required=True)
    p.add_argument("--T", type=float, default=4e5)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--outdir", default="data_paper")
    a = p.parse_args()
    jobs = [dict(src=f, T=a.T, outdir=a.outdir) for f in sorted(glob.glob(a.pattern))]
    print(f"{len(jobs)} continuation job(s), T={a.T:g}")
    with mp.Pool(a.workers) as pool:
        for m in pool.imap_unordered(job, jobs):
            print(m, flush=True)
