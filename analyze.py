"""
Contact-network analysis for the nonreciprocal-colloid modularity test.

Per snapshot:
  Q          Newman modularity of the Louvain partition of the contact network
  nullQ      Q of degree-preserving (double-edge-swap) rewired null, mean of reps
  n_cl       number of connected components with >= 2 particles
  lcf        largest-cluster fraction
Per snapshot pair (consecutive):
  ARI        adjusted Rand index between Louvain partitions (community turnover;
             compare against the same-graph Louvain noise floor, also computed)
  edge_jac   Jaccard distance between contact-edge sets (algorithm-free turnover)
  sigv2      variance of per-particle displacement speeds (order parameter,
             analogue of sigma_v^2 in Hara et al. Fig. 2d/3c)

Usage:
  python analyze.py                 # analyze every .npz in ./data -> results.csv
  python analyze.py --burn 0.5      # discard first half as transient for summary
"""
from __future__ import annotations
import argparse
import glob
import json
import os

import numpy as np
import networkx as nx
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score

CONTACT_FACTOR = 1.1     # edge if r < factor * (s_i + s_j)
NULL_REPS = 3
NULL_EVERY = 4           # compute rewired null every k-th snapshot


def contact_graph(pos, svec, Lbox):
    tree = cKDTree(pos, boxsize=Lbox)
    pairs = tree.query_pairs(CONTACT_FACTOR * 2 * svec.max(), output_type="ndarray")
    G = nx.Graph()
    G.add_nodes_from(range(len(pos)))
    for i, j in pairs:
        d = pos[i] - pos[j]
        d -= Lbox * np.round(d / Lbox)
        if np.hypot(*d) < CONTACT_FACTOR * (svec[i] + svec[j]):
            G.add_edge(int(i), int(j))
    return G


def louvain(G, seed=7):
    comms = nx.community.louvain_communities(G, seed=seed)
    Q = nx.community.modularity(G, comms) if G.number_of_edges() else 0.0
    lab = np.zeros(G.number_of_nodes(), int)
    for k, c in enumerate(comms):
        for n in c:
            lab[n] = k
    return Q, lab


def null_Q(G, reps=NULL_REPS, seed=0):
    qs = []
    for r in range(reps):
        H = G.copy()
        E = H.number_of_edges()
        if E > 10:
            try:
                nx.double_edge_swap(H, nswap=4 * E, max_tries=40 * E, seed=seed + r)
            except nx.NetworkXError:
                pass
        qs.append(louvain(H, seed=seed + r)[0])
    return float(np.mean(qs))


def analyze_file(path, null_every=NULL_EVERY):
    d = np.load(path, allow_pickle=True)
    snaps, svec, Lbox, dts = d["snaps"], d["svec"], float(d["Lbox"]), float(d["dt_snap"])
    T = len(snaps)
    rows = []
    prev_lab, prev_edges = None, None
    for t in range(T):
        G = contact_graph(snaps[t], svec, Lbox)
        Q, lab = louvain(G)
        comps = [c for c in nx.connected_components(G) if len(c) >= 2]
        row = dict(t=t * dts, Q=Q,
                   n_cl=len(comps),
                   lcf=max((len(c) for c in comps), default=0) / len(svec),
                   nullQ=null_Q(G, seed=100 + t) if (t % null_every == 0) else np.nan)
        edges = set(G.edges())
        if prev_lab is not None:
            row["ARI"] = adjusted_rand_score(prev_lab, lab)
            union = edges | prev_edges
            row["edge_jac"] = 1 - len(edges & prev_edges) / len(union) if union else 0.0
            dd = snaps[t] - snaps[t - 1]
            dd -= Lbox * np.round(dd / Lbox)
            row["sigv2"] = float(np.var(np.linalg.norm(dd, axis=1) / dts))
        prev_lab, prev_edges = lab, edges
        rows.append(row)
    # Louvain same-graph noise floor on the final snapshot
    Gf = contact_graph(snaps[-1], svec, Lbox)
    _, l1 = louvain(Gf, seed=7)
    _, l2 = louvain(Gf, seed=99)
    floor = adjusted_rand_score(l1, l2)
    df = pd.DataFrame(rows)
    # chi absent in runs made before the reciprocity knob existed -> those used chi=1
    chi = float(d["chi"]) if "chi" in d.files else 1.0
    meta = dict(file=os.path.basename(path), case=str(d["case"]), seed=int(d["seed"]),
                alpha=float(d["alpha"]), s_ratio=float(d["s_ratio"]), chi=chi,
                ari_floor=floor)
    return df, meta


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--datadir", default="data")
    p.add_argument("--burn", type=float, default=0.5, help="fraction of run discarded as transient")
    p.add_argument("--null-every", type=int, default=NULL_EVERY,
                   help="rewired-null cadence in snapshots; the null is O(E) per rep and "
                        "dominates cost on large graphs")
    p.add_argument("--prefix", default="", help="prefix for output CSVs")
    args = p.parse_args()
    files = sorted(glob.glob(os.path.join(args.datadir, "*.npz")))
    if not files:
        raise SystemExit(f"no .npz files in {args.datadir}; run simulate.py first")
    all_ts, summaries = [], []
    for f in files:
        df, meta = analyze_file(f, args.null_every)
        for k, v in meta.items():
            df[k] = v
        all_ts.append(df)
        late = df.iloc[int(len(df) * args.burn):]
        summaries.append(dict(meta,
                              Q=late.Q.mean(), Q_sd=late.Q.std(),
                              dQ_null=late.Q.mean() - late.nullQ.mean(),
                              ARI=late.ARI.mean(), edge_jac=late.edge_jac.mean(),
                              n_cl=late.n_cl.mean(), lcf=late.lcf.mean(),
                              sigv2=late.sigv2.mean()))
        print(f"{meta['file']}: Q={late.Q.mean():.3f} dQ_null={late.Q.mean()-late.nullQ.mean():.3f} "
              f"ARI={late.ARI.mean():.3f} (floor {meta['ari_floor']:.3f}) "
              f"edge_jac={late.edge_jac.mean():.3f} n_cl={late.n_cl.mean():.1f} "
              f"sigv2={late.sigv2.mean():.2e}", flush=True)
    pd.concat(all_ts).to_csv(f"{args.prefix}timeseries.csv", index=False)
    pd.DataFrame(summaries).to_csv(f"{args.prefix}results.csv", index=False)
    # seed-aggregated summary per condition
    s = pd.DataFrame(summaries)
    agg = s.groupby(["case", "alpha", "s_ratio", "chi"]).agg(
        Q=("Q", "mean"), Q_sem=("Q", "sem"),
        dQ_null=("dQ_null", "mean"),
        ARI=("ARI", "mean"), ARI_sem=("ARI", "sem"),
        edge_jac=("edge_jac", "mean"), edge_jac_sem=("edge_jac", "sem"),
        n_cl=("n_cl", "mean"), sigv2=("sigv2", "mean"), sigv2_sem=("sigv2", "sem"),
        n_seeds=("seed", "count"))
    agg.to_csv(f"{args.prefix}summary_by_condition.csv")
    print("\n== seed-aggregated ==")
    print(agg.to_string())


if __name__ == "__main__":
    main()
