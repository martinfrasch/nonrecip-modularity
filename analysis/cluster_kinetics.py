"""Cluster kinetics from dense-snapshot runs (OPEN_QUESTIONS_PLAN.md, tracks C and D).

For every cluster (size ≥ 2, not the condensate) at snapshot t, find where its members are at
t + 1 and classify:
  split      ≥ 2 pieces of size ≥ 2
  evaporate  one main piece plus monomers only
  merge      the plurality successor also holds ≥ 2 members of another size-≥2 cluster
  absorbed   merge whose successor is the condensate (largest cluster at t+1)
  intact     none of the above
Per-cluster drift ΔS = |successor| − S is reported with absorption excluded (with it included
it reports condensate size, FOLLOWUP.md T2). The condensate itself is tracked for what it
sheds and what it absorbs, size by size — the cycle test.

  python cluster_kinetics.py --pattern "data_paper/*dense500*.npz" "data_dense/*.npz" --out kinetics_events.csv
  python cluster_kinetics.py --summary kinetics_events.csv
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np
import pandas as pd
import networkx as nx
from analysis.analyze import contact_graph

BINS = [2, 3, 5, 7, 10, 14, 20, 30, 50, 100, 10**9]


def labels(pos, sv, L):
    G = contact_graph(pos, sv, L)
    lab = np.full(len(pos), -1); sizes = []
    for k, cc in enumerate(sorted(nx.connected_components(G), key=len, reverse=True)):
        idx = np.fromiter(cc, int); lab[idx] = k; sizes.append(len(idx))
    return lab, np.array(sizes)


def events_for_file(f):
    z = np.load(f)
    sn, sv, L = z["snaps"], z["svec"], float(z["Lbox"])
    dts = float(z["dt_snap"]); chi = float(z["chi"]) if "chi" in z.files else 1.0
    seed = int(z["seed"]); tag = os.path.basename(f)
    lab1, sz1 = labels(sn[0], sv, L)
    rows, cond = [], []
    for t in range(len(sn) - 1):
        lab0, sz0 = lab1, sz1
        lab1, sz1 = labels(sn[t + 1], sv, L)
        # members of each successor cluster, by origin cluster
        for k in range(len(sz0)):
            S = sz0[k]
            if S < 2:
                continue
            mem = np.where(lab0 == k)[0]
            succ, cnt = np.unique(lab1[mem], return_counts=True)
            main = succ[np.argmax(cnt)]
            npieces2 = int(np.sum((cnt >= 2)))
            d = sn[t + 1][mem] - sn[t][mem]; d -= L * np.round(d / L)
            v = np.hypot(*d.mean(0)) / dts
            if k == 0:                              # condensate: what it sheds
                for s_, c_ in zip(succ, cnt):
                    if s_ != main and c_ >= 2:
                        cond.append(dict(file=tag, chi=chi, seed=seed, t=t * dts, kind="shed", size=int(c_)))
                continue
            others = np.where(lab1 == main)[0]
            foreign = others[lab0[others] != k]
            fo_lab, fo_cnt = np.unique(lab0[foreign], return_counts=True) if len(foreign) else ([], [])
            merged = any((c_ >= 2) and (sz0[l_] >= 2) for l_, c_ in zip(fo_lab, fo_cnt))
            absorbed = merged and main == 0
            if npieces2 >= 2:
                kind = "split"
            elif merged:
                kind = "absorbed" if absorbed else "merge"
            elif len(cnt) > 1:
                kind = "evaporate"
            else:
                kind = "intact"
            rows.append(dict(file=tag, chi=chi, seed=seed, t=t * dts, S=int(S), kind=kind,
                             dS=int(sz1[main] - S), v=v, nS=int(np.sum(sv[mem] < 0.14)),
                             keep=int(cnt.max())))
            if absorbed:
                cond.append(dict(file=tag, chi=chi, seed=seed, t=t * dts, kind="absorb", size=int(S)))
    return rows, cond


def summary(ev, cond=None, nboot=500, block=25):
    ev = ev[ev.S >= 2].copy()
    ev["bin"] = pd.cut(ev.S, BINS, right=False)
    out = []
    for (chi, b), g in ev.groupby(["chi", "bin"], observed=True):
        n = len(g)
        if n == 0:
            continue
        na = g.kind != "absorbed"
        bd = g[g.kind.isin(["intact", "evaporate"]) | ((g.kind == "merge") & (g.dS <= 2))]
        out.append(dict(chi=chi, bin=str(b), S_mean=g.S.mean(), n=n, bd_drift=bd.dS.mean(), bd_sem=bd.dS.sem(),
                        p_split=(g.kind == "split").mean(), p_merge=g.kind.isin(["merge", "absorbed"]).mean(),
                        p_absorb=(g.kind == "absorbed").mean(), p_evap=(g.kind == "evaporate").mean(),
                        drift_nonabs=g.dS[na].mean(), drift_nonabs_sem=g.dS[na].std() / np.sqrt(na.sum()),
                        v=g.v.mean()))
    tab = pd.DataFrame(out)

    def crossing(g):
        """Interpolated S (log scale) where the Becker–Döring drift changes sign from − to +, S ≥ 3."""
        g = g[g.S_mean >= 3]
        s, d = np.log(g.S_mean.values), g.bd_drift.values
        for i in range(len(d) - 1):
            if d[i] < 0 <= d[i + 1]:
                return float(np.exp(s[i] + (s[i + 1] - s[i]) * (-d[i]) / (d[i + 1] - d[i])))
        return np.nan

    cross = {}
    rng = np.random.default_rng(0)
    for chi, g in ev.groupby("chi"):
        g = g.copy(); g["blk"] = (g.t // (block * 50)).astype(int).astype(str) + "_" + g.seed.astype(str)
        blks = g.blk.unique()
        base = crossing(summary_bins(g))
        bs = []
        for _ in range(nboot):
            pick = rng.choice(blks, len(blks), replace=True)
            gb = pd.concat([g[g.blk == p] for p in pick])
            bs.append(crossing(summary_bins(gb)))
        bs = np.array(bs); bs = bs[np.isfinite(bs)]
        cross[chi] = (base, np.percentile(bs, 16), np.percentile(bs, 84), len(bs) / nboot)
    return tab, cross


def summary_bins(g):
    """Becker–Döring drift per size bin: monomer/dimer exchange events only (no fission, no fusion)."""
    g = g[g.kind.isin(["intact", "evaporate"]) | ((g.kind == "merge") & (g.dS <= 2))].copy()
    g["bin"] = pd.cut(g.S, BINS, right=False)
    r = g.groupby("bin", observed=True).agg(S_mean=("S", "mean"), bd_drift=("dS", "mean"))
    return r.reset_index()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+")
    p.add_argument("--out", default="results/kinetics_events.csv")
    p.add_argument("--summary", default=None)
    a = p.parse_args()
    if a.summary:
        ev = pd.read_csv(a.summary)
        cond = pd.read_csv(a.summary.replace(".csv", "_condensate.csv"))
        tab, cross = summary(ev)
        pd.set_option("display.width", 200)
        print(tab.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        for chi, (b, lo, hi, frac) in cross.items():
            print(f"chi={chi:g}: Becker–Döring drift zero S* = {b:.2f}  (68% block-bootstrap {lo:.2f}–{hi:.2f}, defined in {frac:.0%} of resamples)")
        for chi, g in cond.groupby("chi"):
            sh, ab = g[g.kind == "shed"].size_, g[g.kind == "absorb"].size_
            print(f"chi={chi:g}: condensate shed {len(sh)} pieces (mean size {sh.mean():.2f}, median {sh.median():.0f}), "
                  f"absorbed {len(ab)} clusters (mean size {ab.mean():.2f}, median {ab.median():.0f})")
        return
    files = sorted(set(f for pat in a.pattern for f in glob.glob(pat)))
    done = set(pd.read_csv(a.out).file) if os.path.exists(a.out) else set()
    for f in files:
        if os.path.basename(f) in done:
            continue
        rows, cond = events_for_file(f)
        pd.DataFrame(rows).to_csv(a.out, mode="a", header=not os.path.exists(a.out), index=False)
        c = a.out.replace(".csv", "_condensate.csv")
        pd.DataFrame(cond).rename(columns={"size": "size_"}).to_csv(c, mode="a", header=not os.path.exists(c), index=False)
        k = pd.Series([r["kind"] for r in rows]).value_counts().to_dict()
        print(f"{os.path.basename(f)[:60]:60s} {len(rows)} cluster-obs  {k}", flush=True)


if __name__ == "__main__":
    main()
