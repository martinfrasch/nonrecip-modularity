"""Thread-parallel force kernel.

The serial kernel walks each pair once (j > i) and writes to BOTH members, which is efficient but
cannot be parallelised without races on F[j]. This version instead gives each cell to one thread
and has it compute the FULL force on its own particles only, looping over all neighbours rather
than half of them. That doubles the pair arithmetic but removes every cross-thread write, so no
atomics or per-thread reduction buffers are needed.

The noise and position update stay serial. They are O(N) and cheap, and keeping them serial
preserves exact reproducibility of the random stream, so a parallel run is comparable with a
serial one at the same seed.
"""
from __future__ import annotations
import numpy as np
from numba import njit, prange
from simulate import RCUT


@njit(cache=True, fastmath=True, parallel=True)
def _forces_par(pos, svec, l2, l4, F, head, nxt, ncell, csize, Lbox, half, alpha, chi):
    """Full force on every particle; cell-parallel, no cross-thread writes."""
    ncell2 = ncell * ncell
    for c in prange(ncell2):
        cx = c // ncell
        cy = c % ncell
        i = head[c]
        while i >= 0:
            fx = 0.0
            fy = 0.0
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    jc = ((cx + dx) % ncell) * ncell + ((cy + dy) % ncell)
                    j = head[jc]
                    while j >= 0:
                        if j != i:
                            rx = pos[i, 0] - pos[j, 0]
                            if rx > half: rx -= Lbox
                            elif rx < -half: rx += Lbox
                            ry = pos[i, 1] - pos[j, 1]
                            if ry > half: ry -= Lbox
                            elif ry < -half: ry += Lbox
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
def _bin(pos, head, nxt, ncell, csize, N):
    head[:] = -1
    for i in range(N):
        c = (int(pos[i, 0] / csize) % ncell) * ncell + (int(pos[i, 1] / csize) % ncell)
        nxt[i] = head[c]
        head[c] = i


def run_par(pos, svec, lvec, nsteps, dt, Lbox, alpha, chi, sigma, seed, snap_every, snaps):
    """Parallel-force integrator; serial noise and update for reproducibility."""
    np.random.seed(seed)
    N = pos.shape[0]
    ncell = int(Lbox / RCUT); csize = Lbox / ncell
    head = np.full(ncell * ncell, -1, np.int64); nxt = np.full(N, -1, np.int64)
    F = np.zeros((N, 2)); l2 = lvec**2; l4 = lvec**4
    half = 0.5 * Lbox
    nstd = sigma * np.sqrt(dt / svec)
    inv = 1.0 / svec
    isnap = 1
    for step in range(1, nsteps + 1):
        _bin(pos, head, nxt, ncell, csize, N)
        _forces_par(pos, svec, l2, l4, F, head, nxt, ncell, csize, Lbox, half, alpha, chi)
        noise = np.random.randn(N, 2)
        pos += (dt * inv)[:, None] * F + nstd[:, None] * noise
        np.mod(pos, Lbox, out=pos)
        if step % snap_every == 0:
            snaps[isnap] = pos
            isnap += 1
    return pos
