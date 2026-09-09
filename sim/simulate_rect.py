"""Rectangular-box integrator, built on the parallel force kernel.

The production kernel assumes a square box: a single Lbox sets the cell sizing, the minimum-image
convention and the periodic wrap. The source experiment's domain is 648 x 360 um, an aspect ratio
of 1.8, and the source-specification replicate approximated it with an equal-area square. This
module removes that approximation so the deviation can be tested rather than assumed.

Everything else -- forces, noise, reproducibility of the random stream -- is identical to
simulate_par.py, so a rectangular run and a square run at matched area, density and seed differ
only in geometry.
"""
from __future__ import annotations
import numpy as np
from numba import njit, prange
from sim.simulate import RCUT


@njit(cache=True, fastmath=True, parallel=True)
def _forces_rect(pos, svec, l2, l4, F, head, nxt, ncx, ncy, csx, csy,
                 Lx, Ly, halfx, halfy, alpha, chi):
    for c in prange(ncx * ncy):
        cx = c // ncy
        cy = c % ncy
        i = head[c]
        while i >= 0:
            fx = 0.0
            fy = 0.0
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    jc = ((cx + dx) % ncx) * ncy + ((cy + dy) % ncy)
                    j = head[jc]
                    while j >= 0:
                        if j != i:
                            rx = pos[i, 0] - pos[j, 0]
                            if rx > halfx: rx -= Lx
                            elif rx < -halfx: rx += Lx
                            ry = pos[i, 1] - pos[j, 1]
                            if ry > halfy: ry -= Ly
                            elif ry < -halfy: ry += Ly
                            r2 = rx * rx + ry * ry
                            if r2 < 1.0:
                                r = np.sqrt(r2)
                                ssum = svec[i] + svec[j]
                                if r < ssum and r > 1e-9:
                                    fm = (ssum - r) / r
                                    fx += fm * rx
                                    fy += fm * ry
                                t = r2 + l2[j]
                                gi = alpha * l4[j] / (t * t * np.sqrt(t))
                                t = r2 + l2[i]
                                gj = alpha * l4[i] / (t * t * np.sqrt(t))
                                gbar = 0.5 * (gi + gj)
                                ci = gbar + chi * (gi - gbar)
                                fx -= ci * rx
                                fy -= ci * ry
                        j = nxt[j]
            F[i, 0] = fx
            F[i, 1] = fy
            i = nxt[i]


@njit(cache=True)
def _bin_rect(pos, head, nxt, ncx, ncy, csx, csy, N):
    head[:] = -1
    for i in range(N):
        c = (int(pos[i, 0] / csx) % ncx) * ncy + (int(pos[i, 1] / csy) % ncy)
        nxt[i] = head[c]
        head[c] = i


def run_rect(pos, svec, lvec, nsteps, dt, Lx, Ly, alpha, chi, sigma, seed, snap_every, snaps):
    np.random.seed(seed)
    N = pos.shape[0]
    ncx = max(3, int(Lx / RCUT)); ncy = max(3, int(Ly / RCUT))
    csx, csy = Lx / ncx, Ly / ncy
    head = np.full(ncx * ncy, -1, np.int64); nxt = np.full(N, -1, np.int64)
    F = np.zeros((N, 2)); l2 = lvec**2; l4 = lvec**4
    halfx, halfy = 0.5 * Lx, 0.5 * Ly
    nstd = sigma * np.sqrt(dt / svec)
    inv = 1.0 / svec
    isnap = 1
    for step in range(1, nsteps + 1):
        _bin_rect(pos, head, nxt, ncx, ncy, csx, csy, N)
        _forces_rect(pos, svec, l2, l4, F, head, nxt, ncx, ncy, csx, csy,
                     Lx, Ly, halfx, halfy, alpha, chi)
        noise = np.random.randn(N, 2)
        pos += (dt * inv)[:, None] * F + nstd[:, None] * noise
        pos[:, 0] = np.mod(pos[:, 0], Lx)
        pos[:, 1] = np.mod(pos[:, 1], Ly)
        if step % snap_every == 0:
            snaps[isnap] = pos
            isnap += 1
    return pos
