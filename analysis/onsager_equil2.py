"""Equilibrate at the Onsager operating point with the STRUCTURALLY DISSIMILAR second channel.

Referee objection: because both nonreciprocal channels were built with the same
pair-co-propulsion structure, the observed antisymmetry of the cross-response may follow from
that construction symmetry rather than from any physical reciprocity relation. The control is to
make the second channel structurally unlike the first -- cubic rather than linear size contrast,
and weighted toward the outer part of the overlap rather than uniform across it (mode2=1 in
onsager.py) -- and ask whether the antisymmetry survives.
"""
from __future__ import annotations
import glob, multiprocessing as mp, os
import numpy as np
from analysis import onsager
from sim.simulate import SIGMA


def equilibrate(a):
    z = np.load(a["src"])
    pos = np.ascontiguousarray(z["snaps"][a["snap"]]).copy()
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    n = int(round(a["T"] / a["dt"]))
    onsager._run2(pos, svec, lvec**2, lvec**4, n, a["dt"], L, float(z["alpha"]),
                  0.5, 0.5, SIGMA, a["seed"], 1)
    out = f"data_onsager2/eq_{a['tag']}.npz"
    np.savez_compressed(out, snaps=pos[None, :, :], svec=svec, lvec=lvec, Lbox=L,
                        alpha=float(z["alpha"]), types=z["types"], dt_snap=1.0, chi=0.5,
                        s_ratio=float(z["s_ratio"]), dt=a["dt"], seed=a["seed"],
                        case=str(z["case"]))
    return out


if __name__ == "__main__":
    os.makedirs("data_onsager2", exist_ok=True)
    srcs = sorted(glob.glob("data/*_x1_*N1000_*seed[12345].npz"))
    starts = [(s, si) for s in srcs for si in np.linspace(60, 99, 8, dtype=int)]
    jobs = [dict(src=s, snap=int(si), T=1e5, dt=0.05, seed=7000 + k, tag=f"long{k:03d}")
            for k, (s, si) in enumerate(starts)]
    print(f"equilibrating {len(jobs)} configurations with the dissimilar channel")
    with mp.Pool(2) as pool:
        for i, _ in enumerate(pool.imap_unordered(equilibrate, jobs)):
            if i % 10 == 0:
                print(f"  {i+1}/{len(jobs)}", flush=True)
    print("done")
