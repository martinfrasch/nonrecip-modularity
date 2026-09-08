"""Configurational estimator of the excess dissipation.

Q̇ = Σ F∘ẋ = Σ_i |F_i|²/s_i + Σ_i D_i ∇_i·F_i   (Stratonovich noise term)
F = F_s + χ F_a, ⟨F_s∘ẋ⟩ = -dU/dt → 0 in stationarity, so
  EPR/N = χ/(T N) [ Σ_i F_a,i·F_i/s_i + Σ_i D_i ∇_i·F_a,i ]
with T = σ²/2, D_i = σ²/(2 s_i). Everything on the right is a function of the configuration.
"""
import sys, glob, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, '/Users/mfrasch/projects/nonrecip-modularity')
from simulate import SIGMA, _run
T = SIGMA**2 / 2

def terms(pos, svec, lvec, L, alpha, chi, eps=0.02):
    N = len(svec)
    pr = cKDTree(pos, boxsize=L).query_pairs(1.0, output_type="ndarray")
    i, j = pr[:, 0], pr[:, 1]
    d = pos[i] - pos[j]; d -= L * np.round(d / L)
    r2 = (d**2).sum(1); r = np.sqrt(r2)
    li, lj = lvec[i], lvec[j]
    ti, tj = r2 + lj**2, r2 + li**2
    gi = alpha * lj**4 / ti**2.5          # acts on i
    gj = alpha * li**4 / tj**2.5          # acts on j
    dgi = -5 * r * alpha * lj**4 / ti**3.5
    dgj = -5 * r * alpha * li**4 / tj**3.5
    gbar = 0.5 * (gi + gj)
    h = -0.5 * (gi - gj)                  # F_a,i = F_a,j = h * d
    dh = -0.5 * (dgi - dgj)
    ssum = svec[i] + svec[j]
    fm = np.where(r < ssum, (ssum - r) / r, 0.0)
    Fs = np.zeros((N, 2)); Fa = np.zeros((N, 2))
    np.add.at(Fs, i, (fm - gbar)[:, None] * d); np.add.at(Fs, j, -(fm - gbar)[:, None] * d)
    np.add.at(Fa, i, h[:, None] * d);            np.add.at(Fa, j, h[:, None] * d)
    F = Fs + chi * Fa
    mu = 1.0 / svec; D = SIGMA**2 / (2 * svec)
    # divergence of F_a: pair (i,j) contributes (2h + r h') to i and -(...) to j
    div = 2 * h + r * dh
    divterm = np.sum(div * (D[i] - D[j]))
    # cutoff delta at r=1: -h(1)(D_i - D_j) per unit r, estimated from a shell
    shell = r > 1 - eps
    delta = -np.sum((h * (D[i] - D[j]))[shell]) / eps
    Qa = np.sum(mu[:, None] * Fa * F) + divterm + delta
    # same for the symmetric part (should be ~0 in stationarity)
    hs = fm - gbar
    dfm = np.where(r < ssum, -ssum / r**2, 0.0)   # d/dr[(ssum-r)/r]
    dhs = dfm - 0.5 * (dgi + dgj)
    divs = 2 * hs + r * dhs
    divs_term = np.sum(divs * (D[i] + D[j]))  # reciprocal pair: same-sign divergence on both
    deltas = -np.sum((hs * (D[i] + D[j]))[shell]) / eps
    Qs = np.sum(mu[:, None] * Fs * F) + divs_term + deltas
    U = np.sum(0.5 * np.where(r < ssum, (ssum - r)**2, 0.0)) + np.sum(  # EHD potential: ∫ gbar r dr
        0.0)  # (EHD potential omitted from U; only its rate matters and Qs carries it)
    return dict(epr=chi * Qa / T / N, Qs=Qs / T / N, Fa2=chi**2 * np.sum(mu[:, None] * Fa * Fa) / T / N,
                cross=chi * np.sum(mu[:, None] * Fa * Fs) / T / N, div=chi * (divterm + delta) / T / N,
                delta=chi * delta / T / N)

def estimate(f, nlast=20, redt=None, reT=50.0):
    z = np.load(f)
    sn, sv, lv, L, a = z["snaps"], z["svec"], z["lvec"], float(z["Lbox"]), float(z["alpha"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    out = []
    idx = np.linspace(len(sn) // 2, len(sn) - 1, nlast, dtype=int)
    for t in idx:
        pos = np.ascontiguousarray(sn[t]).copy()
        if redt is not None:  # re-equilibrate a short window at smaller dt from this snapshot
            ns = int(reT / redt); snaps = np.zeros((2, len(sv), 2)); heat = np.zeros(2)
            _run(pos, sv, lv, lv**2, lv**4, ns, redt, L, a, chi, SIGMA, 777 + t, ns, snaps, heat)
        out.append(terms(pos, sv, lv, L, a, chi))
    keys = out[0].keys()
    return chi, {k: (np.mean([o[k] for o in out]), np.std([o[k] for o in out]) / np.sqrt(len(out))) for k in keys}

if __name__ == "__main__":
    pat = sys.argv[1]; redt = float(sys.argv[2]) if len(sys.argv) > 2 else None
    nlast = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    for f in sorted(glob.glob(pat)):
        chi, r = estimate(f, nlast=nlast, redt=redt)
        print(f"{f.split('/')[-1][:60]:60s} chi={chi:<4g} EPR={r['epr'][0]:+.3e}±{r['epr'][1]:.1e} "
              f"[|Fa|²:{r['Fa2'][0]:+.2e} cross:{r['cross'][0]:+.2e} div:{r['div'][0]:+.2e} delta:{r['delta'][0]:+.1e}] "
              f"Qs={r['Qs'][0]:+.2e}  EPR/chi²={r['epr'][0]/chi**2 if chi else 0:+.3e}", flush=True)
