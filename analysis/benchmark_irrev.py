"""Calibrate the irreversibility estimators against a system with EXACTLY KNOWN entropy
production, at the sample sizes available in the physiological data.

System: biased random walk on a ring of n states (the canonical minimal NESS). Hop right with
probability p, left with q, stay otherwise. The stationary distribution is uniform and the exact
entropy production per step is

    EPR_exact = (p - q) * ln(p / q)

Writing p = r(1+e)/2 and q = r(1-e)/2, this is r*e*ln((1+e)/(1-e)) ~ 2*r*e^2 for small e --
quadratic in the drive e, the same law we measured in the colloid system (EPR = k*chi^2).

An earlier version used a tilted periodic potential; it was discarded because at any barrier
height large enough to be interesting the system is in the activated-crossing regime, where EPR
is exponentially rather than quadratically nonlinear (measured exponent 4.32) and the drift
velocity at low drive is buried in diffusive noise. The ring walk has no such regime.

The question this answers: can a time-series estimator, given only the observable trajectory and
no knowledge of the dynamics, recover the exponent 2 from series as short as one exercise stage?
"""
from __future__ import annotations
import numpy as np
import irreversibility as ir

NSTATE, RATE = 20, 0.5


def simulate(e, nstep, seed=0, obs_noise=0.0):
    rng = np.random.default_rng(seed)
    p, q = RATE * (1 + e) / 2, RATE * (1 - e) / 2
    u = rng.random(nstep)
    step = np.where(u < p, 1, np.where(u < p + q, -1, 0))
    pos = np.cumsum(step) % NSTATE
    x = pos.astype(float)
    if obs_noise:
        x = x + obs_noise * rng.standard_normal(nstep)
    return x


def epr_exact(e):
    p, q = RATE * (1 + e) / 2, RATE * (1 - e) / 2
    return (p - q) * np.log(p / q)


if __name__ == "__main__":
    ES = [0.05, 0.1, 0.2, 0.4]
    print("=== calibration: biased ring walk, EPR_exact = (p-q)ln(p/q) ===")
    le = np.log(ES); lt = np.log([epr_exact(e) for e in ES])
    print(f"  ground-truth exponent of EPR vs drive e : {np.polyfit(le, lt, 1)[0]:.3f}   (theory 2)\n")
    for N in (3000, 300, 150):
        print(f"  --- series length {N} " + ("(whole test)" if N == 3000 else "(one exercise stage)"))
        print("       e     EPR_exact    " + "  ".join(f"{nm+'_z':>11s}" for nm in ir.ESTIMATORS))
        rows = []
        for e in ES:
            x = simulate(e, N, seed=3, obs_noise=0.3)
            zs = {nm: ir.zscore_vs_surrogates(x, nm, nsurr=80, seed=4) for nm in ir.ESTIMATORS}
            rows.append((e, {k: zs[k]["excess"] for k in zs}))
            print(f"     {e:5.2f}   {epr_exact(e):9.5f}    "
                  + "  ".join(f"{zs[nm]['z']:+11.1f}" for nm in ir.ESTIMATORS))
        for nm in ir.ESTIMATORS:
            y = np.array([abs(r[1][nm]) for r in rows])
            if (y > 1e-12).all():
                n = np.polyfit(le, np.log(y), 1)[0]
                flag = "OK" if abs(n - 2) < 0.6 else "biased"
                print(f"       fitted exponent |{nm}| : {n:+.2f}   [{flag}]")
        print()
