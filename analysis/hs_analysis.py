"""Analyse hs_fdt.py twins: integrated FDT ratio and frequency-resolved effective temperature.

Per twin: MSD_x(t) = ⟨Δx_i(t)²⟩ over particles (unperturbed run) and the diagonal response
χ(t) = ⟨ε_i (Δx_i^f(t) − Δx_i^0(t))⟩ / f.  T_eff(t) = MSD(t) / 2χ(t) (= 1 at equilibrium ∀t).
Spectra: C̃(ω) from the velocity autocorrelation C(τ) = MSD''(τ)/2 and R̃(ω) from R(t) = χ'(t),
both by FFT of the finite-difference derivatives with a Hann window; T_eff(ω) = C̃/2R̃'.
Everything pooled over the twelve twins per χ (three seeds × four starts); errors from the
spread over twins.

  python hs_analysis.py --datadir data_hs
"""
from __future__ import annotations
import argparse, glob, os, re
import numpy as np
import pandas as pd
from sim.simulate import SIGMA

T = SIGMA**2 / 2


def load(f):
    z = np.load(f)
    d0, df, eps, ty = z["d0"].astype(float), z["df"].astype(float), z["eps"], z["types"]
    dt = float(z["dt"]) * int(z["rec"]); fval = float(z["f"])
    out = {}
    for lbl, sel in (("L", ty == 0), ("S", ty == 1), ("all", ty >= 0)):
        out[lbl] = dict(msd=(d0[:, sel, 0]**2).mean(1),
                        chi=(eps[sel] * (df[:, sel, 0] - d0[:, sel, 0])).mean(1) / fval)
    return float(z["chi"]), int(z["seed"]), int(z["snap"]), dt, out


def spectra(msd, chi, dt, nmax=None):
    """C(τ) and R(t) on the sampling grid, then one-sided spectra."""
    n = len(msd) if nmax is None else nmax
    msd, chi = msd[:n], chi[:n]
    C = np.gradient(np.gradient(msd, dt), dt) / 2          # velocity autocorrelation, τ ≥ 0
    R = np.gradient(chi, dt)                                # response function, t ≥ 0
    w = np.hanning(2 * n)[n:]                               # half Hann taper on τ ≥ 0
    om = 2 * np.pi * np.fft.rfftfreq(2 * n, dt)
    # C is even in τ: C̃(ω) = 2 ∫_0^∞ C cos ωτ dτ ; R causal: R̃'(ω) = ∫_0^∞ R cos ωt dt
    Ct = 2 * dt * (np.fft.rfft(np.r_[C * w, np.zeros(n)]).real)
    Rt = dt * (np.fft.rfft(np.r_[R * w, np.zeros(n)]).real)
    return om, Ct, Rt


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--datadir", default="data_hs")
    a = p.parse_args()
    files = sorted(glob.glob(os.path.join(a.datadir, "*.npz")))
    per_chi = {}
    for f in files:
        chi, seed, snap, dt, out = load(f)
        per_chi.setdefault(chi, []).append((seed, snap, dt, out))
    rows = []
    tpts = [1, 2, 5, 10, 20, 50, 100, 200, 400, 800]
    print("integrated FDT ratio T_eff(t)/T = MSD/2Tχ, pooled over twins (mean ± sd over twins)")
    for chi in sorted(per_chi):
        tw = per_chi[chi]; dt = tw[0][2]
        print(f"\nchi = {chi}  ({len(tw)} twins)")
        for lbl in ("L", "S", "all"):
            msd = np.mean([o[lbl]["msd"] for _, _, _, o in tw], 0)
            ch = np.mean([o[lbl]["chi"] for _, _, _, o in tw], 0)
            ratio_tw = np.array([o[lbl]["msd"] / (2 * T * o[lbl]["chi"]) for _, _, _, o in tw])
            line = []
            for tp in tpts:
                k = int(round(tp / dt))
                if 0 < k < len(msd):
                    r = msd[k] / (2 * T * ch[k]); sd = np.nanstd(ratio_tw[:, k]) / np.sqrt(len(tw))
                    line.append(f"t={tp:>3}: {r:5.2f}±{sd:4.2f}")
                    rows.append(dict(chi=chi, species=lbl, t=tp, Teff=r, sem=sd, msd=msd[k], chi_resp=ch[k]))
            print(f"  {lbl:3s} " + "  ".join(line))
        # frequency-resolved, pooled over all particles
        msd = np.mean([o["all"]["msd"] for _, _, _, o in tw], 0)
        ch = np.mean([o["all"]["chi"] for _, _, _, o in tw], 0)
        om, Ct, Rt = spectra(msd, ch, dt)
        sel = (om > 0) & (om < np.pi / dt)
        print("  T_eff(ω)/T:", "  ".join(f"ω={om[k]:.3f}: {Ct[k] / (2 * T * Rt[k]):5.2f}"
                                          for k in np.unique(np.geomspace(1, sel.sum() - 1, 9).astype(int))))
        for k in range(1, sel.sum()):
            rows.append(dict(chi=chi, species="all", omega=om[k], Ct=Ct[k], Rt=Rt[k], Teff_w=Ct[k] / (2 * T * Rt[k])))
    pd.DataFrame(rows).to_csv("results/hs_results.csv", index=False)


if __name__ == "__main__":
    main()
