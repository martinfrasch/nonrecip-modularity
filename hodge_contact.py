"""Directed contact graph and Hodge decomposition of the antisymmetric edge flow (track F1).

On every contact edge (r < 1.1(s_i+s_j)) the antisymmetric coefficient w_ij = (a_j − a_i)/2 with
a_i = g(r; l_j), a_j = g(r; l_i) is an antisymmetric edge function: the pair propulsion has
magnitude χ|w|·r and points toward the particle with the larger EHD radius. Because l is
bi-valued, w vanishes on L–L and S–S edges and has one sign on every L–S edge, so the flow is a
two-level potential up to the r-dependence of the kernel. The discrete Hodge decomposition on
the contact graph, w = grad(s) + (cycle-space residual), measures how much of it is *not* a
gradient. Null: the same magnitudes permuted over edges with random signs.

Caveat carried from §4.2 of the paper: this is a Hodge decomposition of a static edge function on
the contact graph, not of the probability current in configuration space.

  python hodge_contact.py --pattern "data_paper/*_N4000_*seed[123].npz" --nsnap 10
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr
from analyze import contact_graph


def edge_flow(pos, sv, lv, L, alpha):
    G = contact_graph(pos, sv, L)
    E = np.array(G.edges(), int)
    i, j = E[:, 0], E[:, 1]
    d = pos[i] - pos[j]; d -= L * np.round(d / L)
    r2 = (d**2).sum(1)
    ai = alpha * lv[j]**4 / (r2 + lv[j]**2)**2.5
    aj = alpha * lv[i]**4 / (r2 + lv[i]**2)**2.5
    return E, 0.5 * (aj - ai), np.sqrt(r2)


def hodge_graph(E, w, N):
    """Least-squares potential s with (grad s)_ij = s_j − s_i; returns cyclic fraction and s."""
    m = len(E)
    B = sp.coo_matrix((np.r_[np.ones(m), -np.ones(m)], (np.r_[np.arange(m), np.arange(m)], np.r_[E[:, 1], E[:, 0]])),
                      shape=(m, N)).tocsr()
    s = lsqr(B, w, atol=1e-12, btol=1e-12, iter_lim=5000)[0]
    res = w - B @ s
    return np.linalg.norm(res) / np.linalg.norm(w), s, res


def one(f, nsnap, rng):
    z = np.load(f)
    sn, sv, lv, L, al = z["snaps"], z["svec"], z["lvec"], float(z["Lbox"]), float(z["alpha"])
    chi = float(z["chi"]) if "chi" in z.files else 1.0
    rows = []
    for t in np.linspace(len(sn) // 2, len(sn) - 1, nsnap, dtype=int):
        E, w, r = edge_flow(sn[t], sv, lv, L, al)
        ls = (lv[E[:, 0]] != lv[E[:, 1]])
        cyc, s, res = hodge_graph(E, w, len(sv))
        # null 1: permute magnitudes over all edges, random sign
        wn = rng.permutation(np.abs(w)) * rng.choice([-1, 1], len(w))
        cyc_null, _, _ = hodge_graph(E, wn, len(sv))
        # null 2: keep the L–S support, permute magnitudes among L–S edges, random sign
        wn2 = np.zeros_like(w); idx = np.where(ls)[0]
        wn2[idx] = rng.permutation(np.abs(w[idx])) * rng.choice([-1, 1], len(idx))
        cyc_null2, _, _ = hodge_graph(E, wn2, len(sv))
        # r-variation: coefficient of variation of |w| on L–S edges
        rows.append(dict(file=os.path.basename(f), chi=chi, t=int(t), n_edges=len(E), frac_LS=ls.mean(),
                         cyclic=cyc, cyclic_null_all=cyc_null, cyclic_null_LS=cyc_null2,
                         cv_w_LS=np.std(np.abs(w[ls])) / np.mean(np.abs(w[ls])),
                         s_gap=(s[lv > 0.14].mean() - s[lv < 0.14].mean())))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--nsnap", type=int, default=10)
    p.add_argument("--out", default="hodge_contact.csv")
    a = p.parse_args()
    rng = np.random.default_rng(1)
    rows = []
    for f in sorted(set(x for pat in a.pattern for x in glob.glob(pat))):
        r = one(f, a.nsnap, rng); rows += r
        d = pd.DataFrame(r)
        print(f"{os.path.basename(f)[:58]:58s} chi={d.chi[0]:<4g} cyclic={d.cyclic.mean():.3f}  "
              f"null(all)={d.cyclic_null_all.mean():.3f} null(LS)={d.cyclic_null_LS.mean():.3f}  cv|w|_LS={d.cv_w_LS.mean():.2f}", flush=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
