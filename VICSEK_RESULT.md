# Chiral Vicsek control: is coarse-grained circulation a signature of biology?

Run 2026-09-02, reproduced 2026-09-05. Code: `vicsek.py`. Reported in `PAPER.md` §3.11.

## Why

§3.11 compares a colloid suspension (no detectable circulation in coarse observables) against
beating cilia (strong circulation). That comparison confounds two contrasts: **biological vs
synthetic**, and **collectively ordered vs spatially distributed and uncoordinated**. A chiral
Vicsek flock is synthetic *and* collectively ordered — an intrinsic turning rate gives it a
collective limit cycle, structurally the same object as the ciliary beat. Applying the identical
signed-area-rate estimator to it separates the two readings.

## Model

    theta_i(t+dt) = arg< e^{i theta_j} >_{j in N_i}  +  omega dt  +  eta * U(-pi, pi)

N=800, L=20, r=1, v=0.3, dt=1, eta=0.25, 6000 steps, first 2000 discarded, 6 seeds per variant.
Nonreciprocity enters as a **vision cone** of half-angle 1.2 rad: i aligns only to neighbours
within that angle of its own heading, so i may see j while j does not see i. Coarse observable is
the mean heading vector (cos theta, sin theta) averaged over particles — matched in spirit to the
two shape modes used for the axoneme. Estimator and phase-randomised surrogate protocol identical
to the colloid and cilia analyses.

## Result

| variant | polar order | median \|z\| | frac \|z\|>2 |
|---|---:|---:|---:|
| standard (reciprocal, no chirality) | 0.709 | 0.55 | 0.00 |
| vision cone (**nonreciprocal**, no chirality) | 0.875 | 0.32 | 0.00 |
| chiral (reciprocal, limit cycle) | 0.015 | 1.80 | 0.33 |
| chiral + vision cone | 0.069 | **5.58** | **1.00** |

Reference points: cilia (single axoneme) median |z| = 3.2, 92% above 2; colloid (single tracked
cluster) median |z| = 0.81–0.94, 0% above 2.

**A purely synthetic system produces circulation stronger than the cilia signal.** Neither
collective order alone nor nonreciprocity alone suffices — both give a null, and the nonreciprocal
non-chiral variant gives a null despite the *highest* polar order of the four. A chiral variant
circulates strongly.

## The signal largely tracks the imposed turning rate

Area rate against omega: 0.0056, 0.0089, 0.0141, 0.0352 for omega = 0.005, 0.01, 0.02, 0.04 —
a ratio to omega of about 0.9 throughout. The measured circulation is therefore largely the
rotation of the mean heading vector at the rate imposed by construction, not an emergent property.

## Consequence for the manuscript

The interpretation that coarse-grained circulation distinguishes biological from non-biological
organisation is **withdrawn**. What the estimator detects is whether the system possesses a
**cyclic collective mode in the chosen projection** — a structural fact about the observable, which
may be emergent (the ciliary beat, from motor coordination) or imposed (a single-particle turning
rate), and which the estimator cannot distinguish.

What survives is narrower and better supported: the colloid suspension has provably positive
entropy production and yet no circulation in any coarse observable, at either the system-averaged
or the single-cluster level. Microscopic irreversibility need not project onto macroscopic
observables. That is a statement about coarse-graining, not about life.

---

*Reproduced 2026-09-05 during the full manuscript audit (`MANUSCRIPT_AUDIT.md`); all four rows
match the values reported in `PAPER.md` §3.11 exactly.*
