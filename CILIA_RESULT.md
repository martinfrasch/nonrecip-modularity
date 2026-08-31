# Does irreversibility survive coarse-graining? Cilia vs the colloid control

First result from the biological arm. Data: [Geyer, Howard & Sartori, Dryad, CC0](https://datadryad.org/dataset/doi:10.5061/dryad.0gb5mkm2j)
— 184 usable *Chlamydomonas* axonemes at 1000 fps across 8 ATP concentrations (50–1000 µM),
272,974 frames. Analysis: `cilia_analysis.py`.

## The test had to be reframed before it could be run

The planned test — is EPR quadratic in the drive — is **ill-posed for this system**. Δμ for ATP
hydrolysis is ≈20 kT at *every* concentration in the series, so cilia are far from equilibrium by
construction, not by measurement. Asking whether they are in linear response answers itself.

The well-posed question is the one our own colloid data set up. That system is nonreciprocal,
dissipative, self-organising and has provably positive entropy production — and yet its
irreversibility **does not survive coarse-graining**: the signed area rate in every coarse
observable plane is indistinguishable from the equilibrium null (all p ≥ 0.38, at two system
sizes). Micro-irreversibility, no macro circulation.

So: **does a biological system differ on exactly that axis?**

## Result

Signed area rate in the plane of the two dominant shape modes (tangent-angle representation,
rigid translation *and rotation* removed), against phase-randomised surrogates:

| [ATP] µM | n | beat Hz | median \|z\| | frac \|z\|>2 | frac \|z\|>3 | \|area/cycle\| |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 10 | 14.5 | 3.3 | 1.00 | 0.70 | 6.07 |
| 66 | 15 | 23.7 | 3.3 | 1.00 | 0.60 | 6.24 |
| 100 | 11 | 27.7 | 3.0 | 0.91 | 0.45 | 6.23 |
| 240 | 19 | 46.3 | 2.7 | 0.89 | 0.32 | 6.15 |
| 370 | 29 | 58.4 | 3.0 | 0.97 | 0.52 | 5.90 |
| 500 | 13 | 91.7 | 4.5 | 0.92 | 0.69 | 5.95 |
| 750 | 44 | 94.8 | 3.2 | 0.89 | 0.55 | 5.71 |
| 1000 | 43 | 65.0 | 3.3 | 0.91 | 0.63 | 6.09 |

**Pooled over 184 axonemes: median |z| = 3.2, 92% exceed |z| = 2, 55% exceed |z| = 3.**
**Colloid control, same estimator on coarse observables: |z| < 1.5, nothing above the null.**

The discriminator works. Both systems are irreversible microscopically; only the biological one
has irreversibility that *survives coarse-graining* into a macroscopic circulation.

## A second, unplanned result

**|area per cycle| is 5.71–6.24 at every ATP level — flat across a 20-fold range of concentration
and a 6-fold range of beat frequency.** That value is ≈2π, the geometric maximum for a circular
orbit in standardised coordinates.

So the beat traces a nearly perfect limit cycle whose *shape* is saturated and ATP-independent.
What ATP sets is the **rate** at which that fixed cycle is traversed: beat frequency rises
14.5 → ~95 Hz, the known Michaelis–Menten behaviour for axonemal dynein.

Circulation rate = (area per cycle, fixed at ≈2π) × (frequency, ATP-dependent). A clean
separation of geometry from kinetics that we did not go looking for.

## Bearing on the Markov-blanket hypothesis

This is direct evidence for Gap 3 of `NOTE_verdict_and_markov_blankets.md`: that the signature of
biological organisation is not directed influence per se — the colloids have plenty — but
**directed influence that persists under coarse-graining**. The colloid is a genuine negative
control for that property, and the cilium is a positive one.

It also reframes the circulation measure. It is a poor probe of nonreciprocity (our strongly
nonreciprocal colloid shows none) but appears to be a good probe of organised, cycle-structured
activity.

## Two bugs found and fixed en route, both silent

- **Rigid rotation was not removed.** A rotated shape has ψ(s) → ψ(s) + const, so the leading
  shape modes captured drift of the axoneme in the field of view rather than the beat. Symptom:
  beat frequency 1–5 Hz, when *Chlamydomonas* axonemes beat at 20–50 Hz. Fixed by subtracting
  the per-frame mean angle before the mean shape.
- **SVD mode signs are arbitrary**, so the sign of any single axoneme's circulation is a
  convention rather than a measurement. Taking a median of signed values across axonemes gave
  randomly alternating results. |z| is the sign-independent statistic.

Both produced plausible-looking output. Neither would have been caught without a physical
expectation (known beat frequency) to check against.

## Caveats

- Circulation in a 2-mode projection is a *lower bound* on irreversibility; the colloid null was
  measured in several planes but neither system's full phase space was searched.
- Phase-randomised surrogates preserve the power spectrum but not higher-order structure; a
  stricter null (IAAFT) would tighten the z-scores.
- Beat frequency is a spectral-peak estimate and is not monotonic at the top of the ATP range
  (91.7, 94.8, then 65.0 Hz) — either estimator noise or genuine heterogeneity, not resolved.
- This shows cilia have macroscopic circulation and the colloid does not. It does **not** yet
  show that this distinction generalises to biological vs non-biological systems at large; two
  systems is an existence proof of a difference, not a law.
