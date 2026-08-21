"""Figures from timeseries.csv / results.csv produced by analyze.py.

  python figures.py                # baseline 4-panel + snapshot render
  python figures.py --sweep       # turnover & Q vs alpha and vs s_II/s_I
"""
from __future__ import annotations
import argparse
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {"monodisperse": "#888888", "bidisperse": "#c0392b"}


def baseline():
    ts = pd.read_csv("timeseries.csv")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    for case, g in ts.groupby("case"):
        c = COLORS.get(case, None)
        m = g.groupby("t").mean(numeric_only=True).reset_index()
        sd = g.groupby("t").std(numeric_only=True).reset_index()
        ax[0, 0].plot(m.t, m.Q, color=c, label=case)
        ax[0, 0].fill_between(m.t, m.Q - sd.Q, m.Q + sd.Q, color=c, alpha=0.2)
        nn = m.dropna(subset=["nullQ"])
        ax[0, 0].plot(nn.t, nn.nullQ, "--", color=c, alpha=0.5, lw=1)
        ax[0, 1].plot(m.t, 1 - m.ARI, color=c, label=case)
        ax[0, 1].fill_between(m.t, 1 - m.ARI - sd.ARI, 1 - m.ARI + sd.ARI, color=c, alpha=0.2)
        fl = 1 - g.ari_floor.mean()
        ax[0, 1].axhline(fl, color=c, ls=":", lw=1)
        ax[1, 0].plot(m.t, m.n_cl, color=c, label=case)
        ax[1, 1].semilogy(m.t, m.sigv2, color=c, label=case)
    ax[0, 0].set(title="Newman modularity Q (dashed: degree-preserving null)", xlabel="t (nondim)", ylabel="Q")
    ax[0, 1].set(title="Community turnover 1\u2212ARI (dotted: same-graph Louvain noise floor)",
                 xlabel="t (nondim)", ylabel="1\u2212ARI")
    ax[1, 0].set(title="Number of clusters (\u22652 particles)", xlabel="t (nondim)", ylabel="n_cl")
    ax[1, 1].set(title="Velocity variance \u03c3\u00b2_v (cf. Hara et al. Fig. 2d)", xlabel="t (nondim)")
    for a in ax.flat:
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
    plt.suptitle("Contact-network modularity test \u2014 Hara et al. PRL 137, 068302 (2026) model, seed-averaged")
    plt.tight_layout()
    plt.savefig("modularity_test.png", dpi=160)
    print("wrote modularity_test.png")

    # final-frame snapshot render, one example per case
    files = {c: sorted(glob.glob(f"data/{c}_*seed1.npz")) or sorted(glob.glob(f"data/{c}_*.npz"))
             for c in ["monodisperse", "bidisperse"]}
    fig2, ax2 = plt.subplots(1, 2, figsize=(11, 5.5))
    for k, case in enumerate(["monodisperse", "bidisperse"]):
        if not files[case]:
            continue
        d = np.load(files[case][0], allow_pickle=True)
        p, s, ty, L = d["snaps"][-1], d["svec"], d["types"], float(d["Lbox"])
        col = np.where(ty == 0, "#2c3e50", "#e67e22")
        for i in range(len(p)):
            ax2[k].add_patch(plt.Circle(p[i], s[i], color=col[i], lw=0))
        ax2[k].set(xlim=(0, L), ylim=(0, L), aspect=1, title=f"{case} (final frame)")
        ax2[k].set_xticks([]); ax2[k].set_yticks([])
    plt.tight_layout()
    plt.savefig("snapshots.png", dpi=160)
    print("wrote snapshots.png")


def sweep():
    s = pd.read_csv("results.csv")
    bi = s[s.case == "bidisperse"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for x, a, ttl in [("alpha", ax[0], "vs nonreciprocal coupling \u03b1\u0302"),
                      ("s_ratio", ax[1], "vs steric ratio s_II/s_I (head\u2192tail-large)")]:
        g = bi.groupby(x).agg(ej=("edge_jac", "mean"), ejs=("edge_jac", "sem"),
                              Q=("Q", "mean"), Qs=("Q", "sem")).reset_index()
        if len(g) < 2:
            a.set_title(f"{ttl} \u2014 no sweep data"); continue
        a.errorbar(g[x], g.ej, yerr=g.ejs, marker="o", color="#c0392b", label="edge turnover")
        a2 = a.twinx()
        a2.errorbar(g[x], g.Q, yerr=g.Qs, marker="s", color="#2c3e50", label="Q")
        a.set(xlabel=x, ylabel="edge Jaccard turnover", title=ttl)
        a2.set_ylabel("Q")
        a.legend(loc="upper left", fontsize=8); a2.legend(loc="lower right", fontsize=8)
        a.grid(alpha=0.3)
    plt.suptitle("Prediction: turnover tracks the antisymmetric coupling; Q stays pinned")
    plt.tight_layout()
    plt.savefig("sweep.png", dpi=160)
    print("wrote sweep.png")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--sweep", action="store_true")
    args = p.parse_args()
    if args.sweep:
        sweep()
    else:
        baseline()
