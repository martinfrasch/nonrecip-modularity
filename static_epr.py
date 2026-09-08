"""Configurational estimator of the excess dissipation (OPEN_QUESTIONS_PLAN.md, track A).

Stratonovich heat  Q̇ = Σ_i F_i∘ẋ_i,  ẋ_i = F_i/s_i + ξ_i,  ⟨g∘ξ_i⟩ = D_i ∇_i·g.
With F = F_s + χ F_a and ⟨F_s∘ẋ⟩ = −dU/dt (≈ 0 in stationarity, measured separately):

    EPR/N = (χ / T N) [ Σ_i F_a,i·F_i / s_i  +  Σ_i D_i ∇_i·F_a,i ],   T = σ²/2, D_i = σ²/(2 s_i)

Everything on the right is a function of the configuration, so the excess dissipation is a
static ensemble average over stored snapshots: no trajectory, no small timestep, no paired
bias probe. The symmetric-sector analogue Q̇_s must vanish in stationarity; on dt = 0.05
snapshots it does not (contact statistics are distorted at O(dt)), so an optional calibration
re-equilibrates a few snapshots at smaller dt and reports the ratio.

  python static_epr.py --pattern "data_paper/*.npz" --out epr_static.csv
  python static_epr.py --pattern "data_paper/*.npz" --out epr_static.csv --calib-dt 0.005 --calib-n 4
"""
from __future__ import annotations
import argparse, glob, multiprocessing as mp, os, re
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from simulate import SIGMA, _run

T = SIGMA**2 / 2


def terms(pos, svec, lvec, L, alpha, chi, eps=0.02):
    """Per-particle rates in units of T per unit time, for one configuration."""
    N = len(svec)
    pr = cKDTree(pos, boxsize=L).query_pairs(1.0, output_type="ndarray")
    i, j = pr[:, 0], pr[:, 1]
    d = pos[i] - pos[j]; d -= L * np.round(d / L)
    r2 = (d**2).sum(1); r = np.sqrt(r2)
    li, lj = lvec[i], lvec[j]
    ti, tj = r2 + lj**2, r2 + li**2
    gi = alpha * lj**4 / ti**2.5                    # acts on i (built from l_j)
    gj = alpha * li**4 / tj**2.5                    # acts on j
    dgi = -5 * r * alpha * lj**4 / ti**3.5
    dgj = -5 * r * alpha * li**4 / tj**3.5
    gbar = 0.5 * (gi + gj)
    h = -0.5 * (gi - gj)                            # F_a,i = F_a,j = h·r_ij  (pair propulsion)
    dh = -0.5 * (dgi - dgj)
    ssum = svec[i] + svec[j]
    fm = np.where(r < ssum, (ssum - r) / r, 0.0)   # reciprocal spring
    hs = fm - gbar                                  # F_s,i = hs·r_ij, F_s,j = −hs·r_ij
    Fs = np.zeros((N, 2)); Fa = np.zeros((N, 2))
    np.add.at(Fs, i, hs[:, None] * d); np.add.at(Fs, j, -hs[:, None] * d)
    np.add.at(Fa, i, h[:, None] * d);  np.add.at(Fa, j, h[:, None] * d)
    F = Fs + chi * Fa
    mu = 1.0 / svec; D = SIGMA**2 / (2 * svec)
    shell = r > 1 - eps                             # cutoff delta at r = 1, shell estimate
    # antisymmetric: ∇_i·(h r_ij) = 2h + r h',  ∇_j·(h r_ij) = −(2h + r h')
    div_a = np.sum((2 * h + r * dh) * (D[i] - D[j])) - np.sum((h * (D[i] - D[j]))[shell]) / eps
    # symmetric (reciprocal): same-sign divergence on both members
    dfm = np.where(r < ssum, -ssum / r**2, 0.0)
    dhs = dfm - 0.5 * (dgi + dgj)
    div_s = np.sum((2 * hs + r * dhs) * (D[i] + D[j])) - np.sum((hs * (D[i] + D[j]))[shell]) / eps
    fa2 = np.sum(mu[:, None] * Fa * Fa); cross = np.sum(mu[:, None] * Fa * Fs)
    Qa = cross + chi * fa2 + div_a
    Qs = np.sum(mu[:, None] * Fs * F) + div_s
    # potential energy (spring + symmetric EHD, zero at the cutoff) for the dU/dt check
    Vl = lambda rr, l: -alpha * l**4 / (3 * (rr**2 + l**2)**1.5)
    U = np.sum(0.5 * np.where(r < ssum, (ssum - r)**2, 0.0)) + np.sum(
        0.5 * ((Vl(r, li) + Vl(r, lj)) - (Vl(1.0, li) + Vl(1.0, lj))))
    return dict(epr=chi * Qa / T / N, a_tilde=fa2 / T / N, b_tilde=(cross + div_a) / T / N,
                Qs=Qs / T / N, U=U / T / N)


def one_file(a):
    f, calib_dt, calib_n, calib_T = a["file"], a["calib_dt"], a["calib_n"], a["calib_T"]
    z = np.load(f)
    sn, sv, lv, L, al = z["snaps"], z["svec"], z["lvec"], float(z["Lbox"]), float(z["alpha"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    dts = float(z["dt_snap"])
    lo = len(sn) // 2
    rows = [terms(sn[t], sv, lv, L, al, chi) for t in range(lo, len(sn))]
    df = pd.DataFrame(rows)
    n = len(df)
    out = dict(file=os.path.basename(f), chi=chi, N=len(sv), seed=int(z["seed"]), L=L,
               dt=float(z["dt"]), nsnap=n, t_window=(n - 1) * dts,
               dUdt=np.polyfit(np.arange(n) * dts, df.U.values, 1)[0])
    for k in ("epr", "a_tilde", "b_tilde", "Qs"):
        out[k] = df[k].mean(); out[k + "_sem"] = df[k].std(ddof=1) / np.sqrt(n)
    if calib_dt:
        idx = np.linspace(lo, len(sn) - 1, calib_n, dtype=int)
        before, after = [], []
        nsamp = a.get("calib_nsamp")            # trajectory mode: many small-dt configurations per start
        for t in idx:
            pos = np.ascontiguousarray(sn[t]).copy()
            before.append(terms(pos, sv, lv, L, al, chi))
            ns = int(calib_T / calib_dt); se = ns // nsamp if nsamp else ns
            snaps = np.zeros((ns // se + 1, len(sv), 2)); heat = np.zeros(len(snaps))
            _run(pos, sv, lv, lv**2, lv**4, ns, calib_dt, L, al, chi, SIGMA, 4242 + int(t), se, snaps, heat)
            first = len(snaps) // 5 if nsamp else len(snaps) - 1     # drop the first fifth as dt-switch transient
            after += [terms(snaps[k], sv, lv, L, al, chi) for k in range(first, len(snaps))]
        calib_n = len(after)
        b, c = pd.DataFrame(before), pd.DataFrame(after)
        out.update(calib_dt=calib_dt, calib_n=calib_n,
                   epr_calib=c.epr.mean(), epr_calib_sem=c.epr.std(ddof=1) / np.sqrt(calib_n),
                   epr_at_calib_snaps=b.epr.mean(), Qs_calib=c.Qs.mean(),
                   a_tilde_calib=c.a_tilde.mean(), b_tilde_calib=c.b_tilde.mean())
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--out", default="epr_static.csv")
    p.add_argument("--calib-dt", type=float, default=None)
    p.add_argument("--calib-n", type=int, default=4)
    p.add_argument("--calib-T", type=float, default=50.0)
    p.add_argument("--calib-nsamp", type=int, default=None,
                   help="trajectory mode: store this many configurations along each small-dt run of length calib-T")
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    files = sorted(set(f for pat in a.pattern for f in glob.glob(pat)))
    done = set(pd.read_csv(a.out).file) if os.path.exists(a.out) else set()
    jobs = [dict(file=f, calib_dt=a.calib_dt, calib_n=a.calib_n, calib_T=a.calib_T, calib_nsamp=a.calib_nsamp)
            for f in files if os.path.basename(f) not in done]
    print(f"{len(files)} files, {len(jobs)} to do, calib_dt={a.calib_dt}", flush=True)
    with mp.Pool(a.workers) as pool:
        for r in pool.imap_unordered(one_file, jobs):
            df = pd.DataFrame([r])
            df.to_csv(a.out, mode="a", header=not os.path.exists(a.out), index=False)
            cal = f"  calib={r['epr_calib']:+.3e} Qs_cal={r['Qs_calib']:+.2f}" if "epr_calib" in r else ""
            print(f"{r['file'][:58]:58s} chi={r['chi']:<4g} N={r['N']:<5d} EPR={r['epr']:+.3e}±{r['epr_sem']:.1e}"
                  f" a~={r['a_tilde']:+.2e} b~={r['b_tilde']:+.2e} Qs={r['Qs']:+.2f}{cal}", flush=True)


if __name__ == "__main__":
    main()
