"""Coarsening law from the structure factor (OPEN_QUESTIONS_PLAN.md, track E).

Domain size L(t) = 2π ∫S(k)dk / ∫k S(k)dk over 0 < k ≤ k_max (Bray's first-moment length),
with S(k) = |Σ_j exp(i k·r_j)|²/N on the periodic k-grid. Runs are chained across
continuations (base → _cont1 → _cont2) so t is absolute. Also records the mass-weighted mean
cluster size ⟨S⟩_w = ΣS²/ΣS and the largest-cluster fraction per snapshot.

  python coarsening_law.py --pattern "data_paper/*.npz" "data_box/*.npz" --out coarsening.csv
"""
from __future__ import annotations
import argparse, glob, os, re
import numpy as np
import pandas as pd
import networkx as nx
from analysis.analyze import contact_graph

KMAX = np.pi          # wavelength ≥ 2 λ-units: coarser than the particle scale


def sk_length(pos, L, nk=48):
    n = np.arange(-nk, nk + 1)
    kx, ky = np.meshgrid(n, n, indexing="ij")
    k = 2 * np.pi / L * np.stack([kx.ravel(), ky.ravel()], 1)
    kmag = np.linalg.norm(k, axis=1)
    keep = (kmag > 0) & (kmag <= KMAX)
    k, kmag = k[keep], kmag[keep]
    # chunk the phase sum to bound memory at N = 22 000
    rho = np.zeros(len(k), complex)
    for s in range(0, len(pos), 2000):
        rho += np.exp(1j * (pos[s:s + 2000] @ k.T)).sum(0)
    S = np.abs(rho)**2 / len(pos)
    return 2 * np.pi * S.sum() / (kmag * S).sum()


def chain(files):
    """Group files by stem, order base < cont1 < cont2 ..., return {stem: [(file, cont_index)]}."""
    groups = {}
    for f in files:
        b = os.path.basename(f)[:-4]
        if "dense" in b:
            continue                              # dense reruns are analysed by their own scripts
        m = re.match(r"(.*?)(?:_cont(\d+))?$", b)
        stem, c = m.group(1), int(m.group(2) or 0)
        groups.setdefault(stem, []).append((f, c))
    return {s: sorted(v, key=lambda x: x[1]) for s, v in groups.items()}


def one_stem(stem, parts, stride):
    rows, t0 = [], 0.0
    for f, c in parts:
        z = np.load(f)
        sn, sv, L = z["snaps"], z["svec"], float(z["Lbox"])
        dts = float(z["dt_snap"]); chi = float(z["chi"]) if "chi" in z.files else 1.0
        start = 1 if c > 0 else 0                 # continuation snapshot 0 duplicates the previous end
        for t in range(start, len(sn), stride):
            G = contact_graph(sn[t], sv, L)
            sizes = np.array([len(cc) for cc in nx.connected_components(G)])
            rows.append(dict(stem=stem, chi=chi, N=len(sv), L=L, seed=int(z["seed"]), cont=c,
                             t=t0 + t * dts, Lk=sk_length(sn[t], L),
                             Sw=(sizes**2).sum() / sizes.sum(), lcf=sizes.max() / len(sv),
                             ncl=int((sizes >= 2).sum())))
        t0 += (len(sn) - 1) * dts
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", nargs="+", required=True)
    p.add_argument("--out", default="results/coarsening.csv")
    p.add_argument("--stride", type=int, default=2)
    a = p.parse_args()
    files = sorted(set(f for pat in a.pattern for f in glob.glob(pat)))
    rows = []
    for stem, parts in chain(files).items():
        r = one_stem(stem, parts, a.stride)
        rows += r
        print(f"{stem[:60]:60s} {len(parts)} part(s) {len(r)} rows  L_k: {r[0]['Lk']:.2f} -> {r[-1]['Lk']:.2f}", flush=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
