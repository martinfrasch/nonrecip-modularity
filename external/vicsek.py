"""Chiral Vicsek control: is coarse-grained circulation a signature of biology, or of
limit-cycle collective organisation?

Section 3.11 compares a colloid suspension (no detectable circulation in coarse observables) with
beating cilia (strong circulation). That comparison confounds two contrasts: biological versus
synthetic, and collectively ordered versus spatially distributed and uncoordinated. A chiral
Vicsek flock is synthetic AND collectively ordered, with an intrinsic turning rate producing a
collective limit cycle -- structurally the same object as the ciliary beat. Applying the identical
signed-area-rate estimator to it discriminates the two readings:

  circulation present  -> the signature is limit-cycle organisation, not biology; Sec 3.11 must
                          be reframed
  circulation absent   -> the biological reading survives a deliberate attempt to break it

Standard Vicsek alignment is reciprocal (mutual metric averaging). We include a vision-cone
variant, in which i aligns only to neighbours within a half-angle phi of its own heading, so that
i may see j while j does not see i -- the natural way to make the alignment nonreciprocal.

Model:
    theta_i(t+dt) = arg< e^{i theta_j} >_{j in N_i}  +  omega dt  +  eta * U(-pi, pi)
    r_i(t+dt)     = r_i(t) + v dt (cos theta_i, sin theta_i)
with omega the intrinsic turning rate (omega = 0 recovers standard Vicsek).
"""
from __future__ import annotations
import argparse
import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def _step(pos, th, L, r, v, dt, eta, omega, cone, noise):
    N = pos.shape[0]
    new = np.empty(N)
    r2 = r * r
    half = 0.5 * L
    for i in range(N):
        sx = np.cos(th[i]); sy = np.sin(th[i])       # include self
        for j in range(N):
            if j == i:
                continue
            dx = pos[j, 0] - pos[i, 0]
            if dx > half: dx -= L
            elif dx < -half: dx += L
            dy = pos[j, 1] - pos[i, 1]
            if dy > half: dy -= L
            elif dy < -half: dy += L
            if dx * dx + dy * dy < r2:
                if cone < 3.15:                       # vision cone -> nonreciprocal
                    d = np.arctan2(dy, dx) - th[i]
                    while d > np.pi: d -= 2 * np.pi
                    while d < -np.pi: d += 2 * np.pi
                    if abs(d) > cone:
                        continue
                sx += np.cos(th[j]); sy += np.sin(th[j])
        new[i] = np.arctan2(sy, sx) + omega * dt + eta * noise[i]
    for i in range(N):
        th[i] = new[i]
        pos[i, 0] = (pos[i, 0] + v * dt * np.cos(th[i])) % L
        pos[i, 1] = (pos[i, 1] + v * dt * np.sin(th[i])) % L


def simulate(N=800, L=20.0, r=1.0, v=0.3, dt=1.0, eta=0.25, omega=0.0,
             cone=np.pi, nsteps=6000, burn=2000, seed=0):
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, (N, 2))
    th = rng.uniform(-np.pi, np.pi, N)
    out = []
    for s in range(nsteps):
        _step(pos, th, L, r, v, dt, eta, omega, cone, rng.uniform(-np.pi, np.pi, N))
        if s >= burn:
            # coarse observables: the two components of the polar order parameter,
            # matched in spirit to the two shape modes used for the axoneme
            out.append((np.cos(th).mean(), np.sin(th).mean()))
    return np.array(out)


def area_rate(x, y, dt=1.0):
    x = (x - x.mean()) / (x.std() + 1e-30)
    y = (y - y.mean()) / (y.std() + 1e-30)
    return 0.5 * np.sum(x[:-1] * y[1:] - y[:-1] * x[1:]) / (len(x) * dt)


def phase_surr(x, rng):
    F = np.fft.rfft(x - x.mean()); ph = rng.uniform(0, 2 * np.pi, len(F)); ph[0] = 0
    if len(x) % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(F) * np.exp(1j * ph), n=len(x)) + x.mean()


def zscore(traj, nsurr=80, seed=0):
    a, b = traj[:, 0], traj[:, 1]
    obs = area_rate(a, b)
    rng = np.random.default_rng(seed)
    null = np.array([area_rate(phase_surr(a, rng), phase_surr(b, rng)) for _ in range(nsurr)])
    sd = null.std(ddof=1)
    return obs, (obs - null.mean()) / sd if sd > 0 else 0.0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=6)
    a = p.parse_args()
    cases = [
        ("standard Vicsek (reciprocal, no chirality)", dict(omega=0.0, cone=np.pi)),
        ("vision cone (nonreciprocal, no chirality)", dict(omega=0.0, cone=1.2)),
        ("chiral (reciprocal, limit cycle)", dict(omega=0.02, cone=np.pi)),
        ("chiral + vision cone (nonreciprocal LC)", dict(omega=0.02, cone=1.2)),
    ]
    print("=== chiral Vicsek circulation control ===")
    print("  cilia (single axoneme): median |z| = 3.2, 92% above |z| = 2")
    print("  colloid (single cluster): median |z| = 0.81-0.94, 0% above |z| = 2\n")
    print("   case                                      polar order   median |z|   frac |z|>2")
    for name, kw in cases:
        zs, ords_ = [], []
        for s in range(a.seeds):
            tr = simulate(seed=s, **kw)
            ords_.append(np.hypot(*tr.mean(0)))
            zs.append(abs(zscore(tr, seed=s)[1]))
        print(f"   {name:41s}  {np.mean(ords_):.3f}      {np.median(zs):7.2f}      {np.mean(np.array(zs) > 2):5.2f}")
