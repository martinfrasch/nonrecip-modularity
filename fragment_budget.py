"""Complete source-sink budget of the fragment population (review round 4, item 5).

Between consecutive dense snapshots, with fragments = clusters of size >= 2 other than the
condensate (largest cluster), count every channel that changes the fragment number:
  +shed        pieces of size >= 2 leaving the condensate
  -absorbed    fragments whose plurality successor is the condensate
  +fission     a fragment splitting into k pieces of size >= 2 adds k-1
  -fusion      a successor fragment fed by m >= 2 fragments removes m-1
  +formation   a successor fragment with no piece of size >= 2 from any single origin cluster
  -dissolution a fragment whose members are all monomers at t+1
and check closure: sum of channels == n_f(t+1) - n_f(t) exactly.

  python fragment_budget.py --pattern "data_dense/*seed1*.npz" "data_paper/*seed1_cont2_dense500.npz"
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np
import pandas as pd
from cluster_kinetics import labels


def budget_file(f):
    z = np.load(f)
    sn, sv, L = z["snaps"], z["svec"], float(z["Lbox"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    lab1, sz1 = labels(sn[0], sv, L)
    rows = []
    for t in range(len(sn) - 1):
        lab0, sz0 = lab1, sz1
        lab1, sz1 = labels(sn[t + 1], sv, L)
        nf0 = int(((sz0 >= 2)).sum() - 1); nf1 = int(((sz1 >= 2)).sum() - 1)
        shed = absorbed = fission = fusion = formation = dissolution = 0
        # origin -> successor member counts
        pairs = pd.DataFrame({"o": lab0, "s": lab1}).value_counts().reset_index(name="n")
        big = pairs[pairs.n >= 2]                              # pieces of size >= 2
        # per origin cluster
        for o, g in pairs.groupby("o"):
            S = sz0[o]
            if S < 2:
                continue
            main = g.loc[g.n.idxmax(), "s"]
            k = int((g.n >= 2).sum())
            if o == 0:                                          # condensate
                shed += int(((g.n >= 2) & (g.s != main)).sum())
                continue
            if main == 0:
                absorbed += 1
            if k >= 2:
                fission += k - 1
            if g.n.max() == 1 and len(g) == S:                  # every member alone
                dissolution += 1
        # per successor cluster: fusion and formation
        for s_, g in big.groupby("s"):
            if s_ == 0 or sz1[s_] < 2:
                continue
            origins = g[g.o.map(lambda o: sz0[o] >= 2)]
            m = origins.o.nunique()
            if m >= 2:
                fusion += m - 1
        for s_ in np.where(sz1 >= 2)[0]:
            if s_ == 0:
                continue
            if not (big.s == s_).any():                          # no piece >= 2 from any origin
                formation += 1
        rows.append(dict(file=os.path.basename(f), chi=chi, t=t, nf0=nf0, nf1=nf1, dnf=nf1 - nf0,
                         shed=shed, absorbed=absorbed, fission=fission, fusion=fusion,
                         formation=formation, dissolution=dissolution))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--out", default="fragment_budget.csv")
    a = p.parse_args()
    rows = []
    for f in sorted(set(x for pat in a.pattern for x in glob.glob(pat))):
        r = budget_file(f); rows += r
        d = pd.DataFrame(r)
        d["sum"] = d.shed - d.absorbed + d.fission - d.fusion + d.formation - d.dissolution
        print(f"{os.path.basename(f)[:58]:58s} chi={d.chi[0]:<4g} closure |sum-dnf| mean {np.abs(d['sum'] - d.dnf).mean():.2f} "
              f"per snapshot: shed {d.shed.mean():.2f} abs {d.absorbed.mean():.2f} fis {d.fission.mean():.2f} "
              f"fus {d.fusion.mean():.2f} form {d.formation.mean():.2f} diss {d.dissolution.mean():.2f}  nf {d.nf0.mean():.1f}", flush=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
