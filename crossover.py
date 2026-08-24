"""Derive the cluster-size crossover where nonreciprocal propulsion stops adding coherently.

A cluster's drive is the vector sum of its pair propulsions p_ij (identical on both members).
If they align, |F| ~ N and speed ~ const. If they are random, |F| ~ sqrt(N) and speed ~ N^-1/2.
The crossover N* is the size of one coherently-driven domain, set by the correlation length
of the propulsion-direction field.
"""
from __future__ import annotations
import glob
import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
from analyze import contact_graph

ALPHA = 0.005


def nonrecip_force(pos, lvec, L, chi=1.0):
    """Per-particle nonreciprocal (propulsive) EHD force; identical on both pair members."""
    F = np.zeros_like(pos)
    pr = cKDTree(pos, boxsize=L).query_pairs(1.0, output_type="ndarray")
    i, j = pr[:, 0], pr[:, 1]
    d = pos[i] - pos[j]; d -= L * np.round(d / L)
    r2 = (d ** 2).sum(1)
    li, lj = lvec[i], lvec[j]
    Gi = ALPHA * lj ** 4 / (r2 + lj ** 2) ** 2.5
    Gj = ALPHA * li ** 4 / (r2 + li ** 2) ** 2.5
    p = (-chi * 0.5 * (Gi - Gj))[:, None] * d
    np.add.at(F, i, p); np.add.at(F, j, p)
    return F


def collect(pattern, chi, stride=4):
    rows = []
    for f in sorted(glob.glob(pattern)):
        z = np.load(f)
        sn, sv, lv, L = z["snaps"], z["svec"], z["lvec"], float(z["Lbox"])
        dts = float(z["dt_snap"])
        for t in range(len(sn) // 2, len(sn) - 1, stride):
            F = nonrecip_force(sn[t], lv, L, chi)
            mag = np.linalg.norm(F, axis=1)
            G = contact_graph(sn[t], sv, L)
            for c in nx.connected_components(G):
                if len(c) < 4:
                    continue
                m = np.fromiter(c, int)
                fnet = float(np.hypot(*F[m].sum(0)))
                fsum = float(mag[m].sum())
                dd = sn[t + 1][m] - sn[t][m]; dd -= L * np.round(dd / L)
                rows.append(dict(N=len(m), fnet=fnet, fsum=fsum,
                                 v=float(np.hypot(*dd.mean(0))) / dts,
                                 zeta=float(sv[m].sum())))
    return rows


def corr_length(pattern, chi, nsnap=3):
    """Spatial correlation of propulsion direction inside the largest cluster."""
    out = []
    for f in sorted(glob.glob(pattern))[:1]:
        z = np.load(f)
        sn, sv, lv, L = z["snaps"], z["svec"], z["lvec"], float(z["Lbox"])
        for t in np.linspace(len(sn) // 2, len(sn) - 2, nsnap, dtype=int):
            F = nonrecip_force(sn[t], lv, L, chi)
            G = contact_graph(sn[t], sv, L)
            big = max(nx.connected_components(G), key=len)
            m = np.fromiter(big, int)
            if len(m) < 100:
                continue
            u = F[m]; n = np.linalg.norm(u, axis=1)
            keep = n > 0
            m, u = m[keep], u[keep] / n[keep, None]
            p = sn[t][m]
            tree = cKDTree(p, boxsize=L)
            for rmax in RBINS:
                pr = tree.query_pairs(rmax, output_type="ndarray")
                if len(pr) < 20:
                    continue
                d = p[pr[:, 0]] - p[pr[:, 1]]; d -= L * np.round(d / L)
                rr = np.hypot(d[:, 0], d[:, 1])
                dot = (u[pr[:, 0]] * u[pr[:, 1]]).sum(1)
                out.append((rmax, rr, dot))
    return out


RBINS = [0.6, 1.2, 2.0, 3.0, 4.5, 6.5]
