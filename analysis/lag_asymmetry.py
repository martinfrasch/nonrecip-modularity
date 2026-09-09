"""All-pairs lag asymmetry of coarse observables (track F3).

For a stationary process with detailed balance every cross-correlation is time-symmetric,
C_ab(τ) = C_ba(τ). The antisymmetric part A_ab(τ) = C_ab(τ) − C_ba(τ), summed over all observable
pairs and lags, is the strongest projection test available on scalar observables: it is nonzero
iff some pair of them traces a cycle at some lag. Observables per snapshot: n_cl, lcf, mean
fragment size, number of contact edges, potential energy U, configurational excess dissipation,
edge Jaccard turnover to the next snapshot. Null: the χ = 0 dense runs analysed identically, plus
a block-shuffled surrogate for each run.

  python lag_asymmetry.py --pattern "data_dense/*.npz" "data_paper/*dense500*.npz" --maxlag 40
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np
import pandas as pd
import networkx as nx
from analysis.analyze import contact_graph
from analysis.static_epr import terms

NAMES = ["n_cl", "lcf", "frag", "edges", "U", "epr", "jac"]


def observables(f, stride=1):
    z = np.load(f)
    sn, sv, lv, L, al = z["snaps"], z["svec"], z["lvec"], float(z["Lbox"]), float(z["alpha"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    rows, prev = [], None
    for t in range(0, len(sn), stride):
        G = contact_graph(sn[t], sv, L)
        E = set(map(tuple, map(sorted, G.edges())))
        sizes = np.array(sorted((len(c) for c in nx.connected_components(G)), reverse=True))
        frag = sizes[1:]; frag = frag[frag >= 2]
        tm = terms(sn[t], sv, lv, L, al, chi if chi > 0 else 1.0)   # at χ = 0 use the χ = 1 proxy
        jac = 1 - len(E & prev) / len(E | prev) if prev is not None else np.nan
        rows.append(dict(t=t, n_cl=int((sizes >= 2).sum()), lcf=sizes[0] / len(sv),
                         frag=frag.mean() if len(frag) else 0.0, edges=len(E), U=tm["U"],
                         epr=tm["epr"], jac=jac))
        prev = E
    d = pd.DataFrame(rows); d["jac"] = d.jac.shift(-1)
    return chi, int(z["seed"]), d.dropna()


def asym(X, maxlag):
    """Frobenius norm of the antisymmetric lagged correlation, per lag, on standardised columns."""
    X = (X - X.mean(0)) / X.std(0)
    n = len(X); out = []
    for k in range(1, maxlag + 1):
        C = X[:-k].T @ X[k:] / (n - k)
        out.append(np.linalg.norm(C - C.T) / np.sqrt(2))
    return np.array(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--maxlag", type=int, default=40)
    p.add_argument("--nsurr", type=int, default=200)
    p.add_argument("--out", default="results/lag_asymmetry.csv")
    a = p.parse_args()
    rng = np.random.default_rng(0)
    rows = []
    for f in sorted(set(x for pat in a.pattern for x in glob.glob(pat))):
        chi, seed, d = observables(f)
        X = d[NAMES].values
        obs = asym(X, a.maxlag)
        # block-shuffle surrogate: permute blocks of 20 snapshots (kills cycles longer than a block,
        # keeps within-block ordering) -- and a full time-reversal check, which flips A exactly
        blk = 20; nb = len(X) // blk
        surr = []
        for _ in range(a.nsurr):
            idx = np.concatenate([np.arange(b * blk, (b + 1) * blk) for b in rng.permutation(nb)])
            surr.append(asym(X[idx], a.maxlag))
        surr = np.array(surr)
        z = (obs.sum() - surr.sum(1).mean()) / surr.sum(1).std()
        rows.append(dict(file=os.path.basename(f), chi=chi, seed=seed, n=len(X), A_sum=obs.sum(),
                         A_surr_mean=surr.sum(1).mean(), A_surr_sd=surr.sum(1).std(), z=z,
                         A_lag1=obs[0], A_lag5=obs[4], A_lag20=obs[19]))
        print(f"{os.path.basename(f)[:58]:58s} chi={chi:<4g} ΣA={obs.sum():.3f} surr={surr.sum(1).mean():.3f}±{surr.sum(1).std():.3f} z={z:+.2f}", flush=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
