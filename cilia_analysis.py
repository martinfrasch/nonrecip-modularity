"""Does irreversibility survive coarse-graining? Cilia vs the colloid control.

The colloid system is nonreciprocal, dissipative and self-organising, with provably positive
entropy production -- yet its irreversibility does NOT survive coarse-graining: the signed area
rate in every coarse observable plane is indistinguishable from the equilibrium null (all
p >= 0.38, at two system sizes). Micro-irreversibility, no macro circulation.

This asks the same question of a biological system. Cilia beat as a limit cycle, so if
biological organisation is characterised by irreversibility that survives coarse-graining
(the Gap-3 hypothesis in NOTE_verdict_and_markov_blankets.md), the coarse-grained circulation
should be large where the colloid's is zero.

Note what is NOT being tested: whether EPR is quadratic in the ATP drive. Delta_mu for ATP
hydrolysis is ~20 kT at every concentration here, so cilia are far from equilibrium by
construction and the linear-response question is ill-posed for them.

Shape representation: tangent angle psi(s) along the arclength, which removes rigid-body
translation and rotation, then PCA over frames -> two dominant shape modes.
"""
from __future__ import annotations
import numpy as np, scipy.io as sio, glob, re, os


def tangent_angles(frames, npts=50):
    """(x,y) contours -> unwrapped tangent angle psi(s), resampled to npts, per frame."""
    out = []
    for fr in frames:
        arr = np.asarray(fr, float)
        if arr.ndim != 2 or arr.shape[0] < 10 or arr.shape[1] < 2:
            out.append(None); continue      # untracked / malformed frame
        xy = arr[:, :2]
        xy = xy[np.isfinite(xy).all(1)]
        if len(xy) < 10:
            out.append(None); continue
        d = np.diff(xy, axis=0)
        s = np.concatenate([[0], np.cumsum(np.hypot(d[:, 0], d[:, 1]))])
        psi = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
        sm = 0.5 * (s[:-1] + s[1:])
        out.append(np.interp(np.linspace(sm[0], sm[-1], npts), sm, psi))
    good = [o for o in out if o is not None]
    return np.array(good) if good else np.empty((0, npts))


def area_rate(x, y, dt):
    """Signed area swept per unit time in a standardised 2-D observable plane.
    Exactly zero under detailed balance for any observable pair."""
    x = (x - x.mean()) / (x.std() + 1e-30)
    y = (y - y.mean()) / (y.std() + 1e-30)
    return 0.5 * np.sum(x[:-1] * y[1:] - y[:-1] * x[1:]) / (len(x) * dt)


def phase_surrogate(x, rng):
    F = np.fft.rfft(x - x.mean()); ph = rng.uniform(0, 2*np.pi, len(F)); ph[0] = 0
    if len(x) % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(F)*np.exp(1j*ph), n=len(x)) + x.mean()


def iaaft_surrogate(x, rng, niter=30):
    """Iterative amplitude-adjusted Fourier transform (Schreiber & Schmitz 1996): preserves both the
    power spectrum and the amplitude distribution of x. A stricter null than phase randomisation."""
    xs = np.sort(x); amp = np.abs(np.fft.rfft(x)); y = rng.permutation(x)
    for _ in range(niter):
        y = np.fft.irfft(amp * np.exp(1j * np.angle(np.fft.rfft(y))), n=len(x))
        y = xs[np.argsort(np.argsort(y))]
    return y


SURROGATE = phase_surrogate


def analyse_axoneme(ax, nsurr=30, seed=0):
    psi = tangent_angles(np.atleast_1d(ax.XY_Data))
    if len(psi) < 200:
        return None
    dt = float(ax.dt)
    # remove rigid-body ROTATION: a rotated shape has psi(s) -> psi(s) + const, so the
    # per-frame mean angle is pure rotation and must go before the mean shape is subtracted.
    # Without this the leading modes capture drift of the axoneme in the field of view
    # rather than the beat, and the beat frequency is badly underestimated.
    psi = psi - psi.mean(1, keepdims=True)
    psi = psi - psi.mean(0)
    u, s, vt = np.linalg.svd(psi, full_matrices=False)
    a1, a2 = u[:, 0]*s[0], u[:, 1]*s[1]
    var2 = (s[0]**2 + s[1]**2) / (s**2).sum()
    obs = area_rate(a1, a2, dt)
    rng = np.random.default_rng(seed)
    null = np.array([area_rate(SURROGATE(a1, rng), SURROGATE(a2, rng), dt)
                     for _ in range(nsurr)])
    sd = null.std(ddof=1)
    # beat frequency from the dominant spectral peak of mode 1
    f = np.fft.rfftfreq(len(a1), dt); P = np.abs(np.fft.rfft(a1 - a1.mean()))**2
    freq = f[1:][np.argmax(P[1:])]
    return dict(n=len(psi), var2=var2, area=obs, null=null.mean(), null_sd=sd,
                z=(obs - null.mean())/sd if sd > 0 else 0.0, freq=freq,
                area_per_cycle=obs/freq if freq > 0 else np.nan)


if __name__ == "__main__":
    import sys
    if "--iaaft" in sys.argv:
        SURROGATE = iaaft_surrogate
        print("  surrogate: iAAFT (spectrum + amplitude distribution preserved)")
    files = sorted(glob.glob('data_cilia/**/WT_*uM-ATP.mat', recursive=True),
                   key=lambda p: int(re.search(r'WT_(\d+)uM', p).group(1)))
    print("=== coarse-grained circulation in cilia shape space ===")
    print("  colloid control: area-rate z indistinguishable from zero in every plane (all p >= 0.38)\n")
    print("  NOTE: SVD mode signs are arbitrary, so the SIGN of a single axoneme's circulation")
    print("  is a convention, not a measurement. |z| is the sign-independent statistic.\n")
    ALL = []
    print("  [ATP]uM  n_ax  beat Hz   median|z|   frac |z|>2   frac |z|>3   |area/cycle|")
    for path in files:
        atp = int(re.search(r'WT_(\d+)uM', path).group(1))
        dp = np.atleast_1d(sio.loadmat(path, squeeze_me=True, struct_as_record=False)['datapool'])
        R = [r for r in (analyse_axoneme(a, seed=i) for i, a in enumerate(dp)) if r]
        if not R: continue
        z = np.abs(np.array([r['z'] for r in R]))
        apc = np.abs(np.array([r['area_per_cycle'] for r in R]))
        print(f"   {atp:6d}  {len(R):4d}   {np.median([r['freq'] for r in R]):6.1f}    "
              f"{np.median(z):7.1f}     {np.mean(z>2):7.2f}      {np.mean(z>3):7.2f}      {np.median(apc):7.2f}")
        ALL.extend(z)
    A = np.array(ALL)
    print(f"\n  POOLED over all {len(A)} axonemes: median |z| = {np.median(A):.1f}, "
          f"{np.mean(A>2)*100:.0f}% exceed |z|=2, {np.mean(A>3)*100:.0f}% exceed |z|=3")
    print(f"  COLLOID CONTROL (same estimator, coarse observables): |z| < 1.5, all p >= 0.38")
