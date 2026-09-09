"""Frequency-resolved effective temperature on single-particle coordinates (track B).

Apply a random-sign force f ε_i x̂ (ε_i = ±1 i.i.d.) to every particle and run with the same
noise stream as an unperturbed twin. The paired response of Σ_i ε_i x_i / (f N) is the mean
DIAGONAL displacement response χ(t) (cross responses average out over ε), which is the quantity
conjugate to the mean single-particle MSD. At equilibrium MSD(t) = 2T χ(t) at every t; out of
equilibrium T_eff(t) = MSD(t)/2χ(t), and the Fourier transforms of the velocity autocorrelation
C(τ) = MSD''/2 and of R(t) = χ'(t) give T_eff(ω) = C̃(ω)/2R̃'(ω) — the Cugliandolo–Kurchan object.
The Harada–Sasa violation spectrum V(ω) = [C̃(ω) − 2T R̃'(ω)]/μ, integrated over ω, is the
per-particle dissipation, to be compared with the configurational estimate of static_epr.py.

  python hs_fdt.py --pattern "data_paper/*_x1_N4000*seed1.npz" --T 1000 --f 1e-4 --rec 10 --outdir data_hs
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp, os
import numpy as np
from numba import njit
from sim.simulate import SIGMA, RCUT


@njit(cache=True, fastmath=True)
def _run_rec(pos, svec, l2, l4, fext, nsteps, dt, Lbox, alpha, chi, sigma, seed, rec, out):
    """Integrate; store unwrapped displacement of every particle every `rec` steps in out[k]."""
    np.random.seed(seed)
    N = pos.shape[0]
    ncell = int(Lbox / RCUT); csize = Lbox / ncell
    head = np.full(ncell * ncell, -1, np.int64); nxt = np.full(N, -1, np.int64)
    F = np.zeros((N, 2)); disp = np.zeros((N, 2))
    half = 0.5 * Lbox
    nstd = np.empty(N)
    for i in range(N):
        nstd[i] = sigma * np.sqrt(dt / svec[i])
    k = 1
    for step in range(1, nsteps + 1):
        head[:] = -1
        for i in range(N):
            c = (int(pos[i, 0] / csize) % ncell) * ncell + (int(pos[i, 1] / csize) % ncell)
            nxt[i] = head[c]; head[c] = i
        for i in range(N):
            F[i, 0] = fext[i]; F[i, 1] = 0.0
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
                                    if rx > half: rx -= Lbox
                                    elif rx < -half: rx += Lbox
                                    ry = pos[i, 1] - pos[j, 1]
                                    if ry > half: ry -= Lbox
                                    elif ry < -half: ry += Lbox
                                    r2 = rx * rx + ry * ry
                                    if r2 < 1.0:
                                        r = np.sqrt(r2); ssum = svec[i] + svec[j]
                                        if r < ssum and r > 1e-9:
                                            fm = (ssum - r) / r
                                            F[i, 0] += fm * rx; F[i, 1] += fm * ry
                                            F[j, 0] -= fm * rx; F[j, 1] -= fm * ry
                                        t = r2 + l2[j]
                                        gi = alpha * l4[j] / (t * t * np.sqrt(t))
                                        t = r2 + l2[i]
                                        gj = alpha * l4[i] / (t * t * np.sqrt(t))
                                        gbar = 0.5 * (gi + gj)
                                        ci = gbar + chi * (gi - gbar)
                                        cj = gbar + chi * (gj - gbar)
                                        F[i, 0] -= ci * rx; F[i, 1] -= ci * ry
                                        F[j, 0] += cj * rx; F[j, 1] += cj * ry
                                j = nxt[j]
                    i = nxt[i]
        for i in range(N):
            inv = 1.0 / svec[i]
            ddx = dt * inv * F[i, 0] + nstd[i] * np.random.randn()
            ddy = dt * inv * F[i, 1] + nstd[i] * np.random.randn()
            disp[i, 0] += ddx; disp[i, 1] += ddy
            pos[i, 0] = (pos[i, 0] + ddx) % Lbox
            pos[i, 1] = (pos[i, 1] + ddy) % Lbox
        if step % rec == 0:
            for i in range(N):
                out[k, i, 0] = disp[i, 0]; out[k, i, 1] = disp[i, 1]
            k += 1


def job(a):
    z = np.load(a["src"])
    sv, lv, L, al = z["svec"], z["lvec"], float(z["Lbox"]), float(z["alpha"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    seed = int(z["seed"]); dt = float(z["dt"])
    tag = f"{os.path.basename(a['src'])[:-4]}_snap{a['snap']}_f{a['f']:g}"
    out = os.path.join(a["outdir"], tag + ".npz")
    if os.path.exists(out):
        return "skip " + tag
    rng = np.random.default_rng(seed * 31 + a["snap"])
    eps = rng.choice([-1.0, 1.0], len(sv))
    ns = int(a["T"] / dt); nk = ns // a["rec"] + 1
    res = {}
    for lbl, f in (("0", 0.0), ("f", a["f"])):
        pos = np.ascontiguousarray(z["snaps"][a["snap"]]).copy()
        d = np.zeros((nk, len(sv), 2))
        _run_rec(pos, sv, lv**2, lv**4, f * eps, ns, dt, L, al, chi, SIGMA, 555 + seed, a["rec"], d)
        res[lbl] = d.astype(np.float32)
    np.savez_compressed(out, d0=res["0"], df=res["f"], eps=eps, svec=sv, lvec=lv, types=z["types"],
                        chi=chi, seed=seed, f=a["f"], dt=dt, rec=a["rec"], T=a["T"], snap=a["snap"])
    return "done " + tag


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--T", type=float, default=1000.0)
    p.add_argument("--f", type=float, default=1e-4)
    p.add_argument("--rec", type=int, default=10)
    p.add_argument("--snaps", default="-1", help="comma list of snapshot indices to start from")
    p.add_argument("--outdir", default="data_hs")
    p.add_argument("--workers", type=int, default=3)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    srcs = sorted(set(x for pat in a.pattern for x in glob.glob(pat)))
    jobs = [dict(src=s, snap=int(k), T=a.T, f=a.f, rec=a.rec, outdir=a.outdir)
            for s in srcs for k in a.snaps.split(",")]
    print(f"{len(jobs)} jobs", flush=True)
    with mp.Pool(a.workers) as pool:
        for m in pool.imap_unordered(job, jobs):
            print(m, flush=True)
