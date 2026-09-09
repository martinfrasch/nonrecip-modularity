"""Equilibrate configurations AT the Onsager operating point before measuring the response.

The first attempt measured short windows starting from configurations equilibrated at
(chi=1, chi2=0). Those windows therefore captured relaxation toward the (0.5, 0.5) steady
state rather than the steady-state response, and L21 drifted monotonically with window length
(-1.76e-5 at T=50 to -3.84e-5 at T=400, no plateau). The apparent reciprocity violation was an
artifact of that transient.

This equilibrates at the operating point first, then re-measures.
"""
from __future__ import annotations
import glob, multiprocessing as mp, os
import numpy as np
from analysis.onsager import _run2
from sim.simulate import SIGMA


def equilibrate(a):
    z = np.load(a["src"])
    pos = np.ascontiguousarray(z["snaps"][a["snap"]]).copy()
    svec, lvec, L = z["svec"], z["lvec"], float(z["Lbox"])
    n = int(round(a["T"] / a["dt"]))
    _run2(pos, svec, lvec**2, lvec**4, n, a["dt"], L, float(z["alpha"]),
          a["chi"], a["chi2"], SIGMA, a["seed"])
    out = f"data_onsager/eq_{a['tag']}.npz"
    np.savez_compressed(out, snaps=pos[None, :, :], svec=svec, lvec=lvec,
                        Lbox=L, alpha=float(z["alpha"]), types=z["types"],
                        dt_snap=1.0, chi=a["chi"], s_ratio=float(z["s_ratio"]),
                        dt=a["dt"], seed=a["seed"], case=str(z["case"]))
    return out


if __name__ == "__main__":
    os.makedirs("data_onsager", exist_ok=True)
    srcs = sorted(glob.glob("data/*_x1_*N1000_*seed[12345].npz"))
    jobs = [dict(src=s, snap=int(si), chi=0.5, chi2=0.5, T=2e4, dt=0.05,
                 seed=3000 + k, tag=f"{k:03d}")
            for k, (s, si) in enumerate([(s, si) for s in srcs
                                         for si in np.linspace(60, 99, 8, dtype=int)])]
    print(f"equilibrating {len(jobs)} configurations at (chi, chi2) = (0.5, 0.5), T=2e4")
    with mp.Pool(8) as pool:
        for i, o in enumerate(pool.imap_unordered(equilibrate, jobs)):
            if i % 10 == 0: print(f"  {i+1}/{len(jobs)}", flush=True)
    print("done")
