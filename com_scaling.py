"""Does the ballistic centre-of-mass drift survive the thermodynamic limit?

Sec 3.5 shows pair propulsions add incoherently, C ~ N^-0.53. Applied to the whole system that
gives net internal force ~ sqrt(N), hence V_com ~ N^-1/2 and <dR^2_com> ~ t^2/N: the exponent
would be 2 at every finite N while the effect vanishes as N grows. Reviewer objection; testable
from the existing trajectories.

Also: if the drift DIRECTION decorrelates, t^~2 is window-limited and the coordinate is
diffusive with large D rather than genuinely ballistic.
"""
import glob, os
import numpy as np

def com_track(path):
    d = np.load(path, allow_pickle=True)
    sn = d["snaps"]; L = float(d["Lbox"]); dts = float(d["dt_snap"])
    chi = float(d["chi"]) if "chi" in d.files else 1.0
    N = len(d["svec"])
    dd = np.diff(sn, axis=0)                    # (T-1, N, 2)
    dd -= L * np.round(dd / L)                  # minimum image
    step = dd.mean(axis=1)                      # COM displacement per snapshot
    return np.cumsum(step, axis=0), dts, N, chi, L

def msd_exponent(traj, dts, nlag=8):
    T = len(traj)
    lags = np.unique(np.logspace(0, np.log10(T // 3), nlag).astype(int))
    m = []
    for lg in lags:
        d = traj[lg:] - traj[:-lg]
        m.append((d**2).sum(1).mean())
    m = np.array(m); lags = np.array(lags, float)
    return np.polyfit(np.log(lags * dts), np.log(m), 1)[0], lags * dts, m

print(f"{'file':52s} {'N':>6} {'chi':>4} {'|V_com|':>10} {'V*sqrt(N)':>10} {'MSD exp':>8} {'dir autocorr':>12}")
rows = []
for pat in ["data_paper/bidisperse_*_x1_N4000_L24_T400000_seed?.npz",
            "data_paper/bidisperse_*_x1.5_N4000_L24_T400000_seed?.npz",
            "data_box/*_x1_N9000_*cont1.npz", "data_box/*_x1.5_N9000_*cont1.npz",
            "data_box/*_x1.5_N16000_*cont1.npz",
            "data_validate/*_cont1.npz"]:
    for f in sorted(glob.glob(pat)):
        if "cont" in f and "data_paper" in f: continue
        tr, dts, N, chi, L = com_track(f)
        tot_t = len(tr) * dts
        V = np.linalg.norm(tr[-1]) / tot_t
        ex, lg, m = msd_exponent(tr, dts)
        # direction autocorrelation: unit velocity over halves of the run
        h = len(tr) // 2
        v1 = tr[h] - tr[0]; v2 = tr[-1] - tr[h]
        ac = float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-30))
        print(f"{os.path.basename(f)[:52]:52s} {N:6d} {chi:4g} {V:10.3e} {V*np.sqrt(N):10.3e} {ex:8.3f} {ac:12.3f}")
        rows.append((N, chi, V, ex, ac))

import collections
print("\n=== aggregated by (N, chi) ===")
g = collections.defaultdict(list)
for N, chi, V, ex, ac in rows: g[(N, chi)].append((V, ex, ac))
for k in sorted(g):
    a = np.array(g[k])
    print(f"  N={k[0]:6d} chi={k[1]:<4g} n={len(a)}  |V_com|={a[:,0].mean():.3e}  "
          f"V*sqrt(N)={a[:,0].mean()*np.sqrt(k[0]):.3e}  MSDexp={a[:,1].mean():.2f}  dir_ac={a[:,2].mean():+.2f}")
print("\n=== V_com vs N at fixed chi ===")
for chi in (1.0, 1.5):
    pts = [(k[0], np.mean([x[0] for x in g[k]])) for k in sorted(g) if k[1] == chi]
    if len(pts) >= 2:
        Ns = np.array([p[0] for p in pts], float); Vs = np.array([p[1] for p in pts])
        sl = np.polyfit(np.log(Ns), np.log(Vs), 1)[0]
        print(f"  chi={chi}: N={list(Ns.astype(int))}  V={[f'{v:.2e}' for v in Vs]}   slope = {sl:+.3f}  (expect -0.5 if incoherent)")
