# Note: the variational family, the decompression vision, and a biological tier test

Saved 2026-08-30 for later work. Companion to `TODO.md` (the ordered programme) and
`VARIATIONAL_FAMILY.md` (the technical placement).

---

## 1. Yes — Onsager holds here, firmly

EPR = kχ² with k ≈ 2.93×10⁻³, **flat to 4.8%** across χ = 1.5–8: a fourfold range of drive, a
fifteenfold range of EPR, four independent chances to fail. The system is in the linear-response
tier across its entire accessible range, and a Rayleighian is the right variational object for
its dynamics.

## 2. The negative result is the more useful one

This system has nonreciprocity, sustained dissipation, spontaneous structure formation,
long-lived organised states, and a nucleation-separated two-phase architecture — every ingredient
usually invoked as a hallmark of active matter. **And it is still Onsager.**

So nonreciprocity + dissipation + self-organisation is *not sufficient* to leave linear response.
Whatever distinguishes living organisation, it is not merely being driven and non-reciprocal.

That makes this system a **control for "driven but not alive"**: any claimed biological signature
must distinguish itself from *this*, not merely from equilibrium. A much harder and more useful
bar than the usual one.

## 3. Where the Triple-Action sits

NWAP's I_max − E_min + symmetry has the structure of **Maximum Caliber** (Jaynes) — maximise path
entropy subject to dynamical constraints. Max Caliber contains an information term by
construction, reduces to MaxEnt statically and to Onsager near equilibrium, and is the leading
candidate principle for far-from-equilibrium dynamics.

That placement would make NWAP a **constrained member of a principled family** — constraints
network-structural rather than thermodynamic — instead of a bespoke functional. It also predicts
checkable limits: NWAP should reduce to a Rayleighian at weak drive and to free-energy
minimisation at zero drive. Either failure is informative.

**What keeps it falsifiable:** each tier carries a signature and membership must be *measured*.
A1 quadratic dissipation (run, passed), A2 Onsager reciprocity (not run, needs a second drive),
A3 single effective temperature (not run). A framework claiming "some variational principle
applies" forbids nothing; one claiming tier II makes three predictions that can each fail.

**The boundary our own data mark:** EPR is quadratic *at fixed structure*, but structural
response is strongly non-linear (cluster count dips at χ=0.25, then rises tenfold). The
Rayleighian governs dynamics at fixed structure; **structure selection is uncovered by any
current variational principle.** That gap is where new theory has to go.

## 4. Horizontal scales: the method this project produced

Measure the **rates** of upward flow (merge/growth) and downward flow (split/fragmentation)
versus the scale variable, then classify crossings:

> **An attractor in the scale-flow is a populated level of organisation.
> A repeller is a boundary between levels.**

Demonstrated: split/merge rates cross at S\* ≈ 7–8 particles, and the crossing is a **repeller** —
creating two levels while populating neither at S\* itself.

**Caution learned three times:** a coarse-grained increment is not a rate. Net drift over the
production interval returned +280 to +1284 in *every* size bin — silent nonsense with plausible
magnitudes and the wrong sign structure.

## 5. The decompression vision has a rigorous home

The metaphysical framing is not testable. The structure inside it is, and it already has a name.

**A scale-to-scale map with fixed points classified by stability is a renormalisation group
flow.** Decompression from a compressed specification is an *inverse* RG flow, and the horizontal
scales are its fixed points — attractors where organisation accumulates, repellers where one
level ends and the next begins. Our split/merge measurement is an empirical RG fixed-point
analysis in miniature.

**The falsifiable core:** if organisation at scale n+1 is decompressed from scale n, then n+1 must
be largely predictable from n by a scale-local rule. Measurable as conditional description length
across adjacent scales — H(n+1 | n) ≪ H(n+1) indicates decompression-like generation. Estimable
by compression-based methods, which is fitting given where the intuition came from.

**The null model is what makes it science:** the same measurement on a system with matched
marginal statistics but scrambled cross-scale structure. Without it the quantity always looks
impressive. *Keep the vision; publish the number.*

## 6. The vagus postulate is the strongest concrete piece

If the vagus is a vertical organiser, cross-scale directed information should be concentrated
through it — so transfer entropy between scale levels should collapse under vagal blockade while
within-level dynamics persist. A paired-recording experiment with a clear negative outcome.

**The two threads converge:** the reciprocity test (A2) that would settle the tier question in
physics needs a second independent drive, and vagal blockade *is* a second drive in physiology.
Worth designing as one experiment, not two.

## 7. The biological tier test — data now in hand

**The answerable question:** not "does a variational principle exist for living systems" but
**"at which tier do they sit, and does anything biological fail the tier-II tests this colloid
passes?"**

**Prediction:** biological systems show *non-quadratic* EPR scaling in their drive. If a driven
biological system came out quadratic, it would be in the same class as this colloid suspension —
organised, dissipative, and unremarkable.

### Dataset: downloaded and verified

[ACTES graded cycloergometer exercise testing, PhysioNet v1.0.0](https://physionet.org/content/actes-cycloergometer-exercise/1.0.0/)
→ `data_actes/`, all four files SHA-256 verified against the published manifest.

| | |
|---|---|
| subjects | 18 adolescent athletes (fencing 10, kayak 6, triathlon 2), ages 12–18 |
| beats | 52,062 total; 2,000–3,854 per subject (median 2,750) |
| duration | 947–1,673 s per test |
| columns | ID, time (s), **RR (ms)**, VO2 (L/min), **power (W)** |
| peak power | 140–335 W |
| per subject | age, weight, height, **P_vt1, P_vt2** (ventilatory thresholds), sport |

**Verified suitable:** power is genuinely **stepwise graded** — subject 1 shows 10 levels with
≥20 beats (0 W for 987 s baseline, then 50, 65, 80, 95, 110, 125, 140, 155 W in ~60 s steps).
That is a clean drive ladder, directly parallel to our χ sweep.

**The sample-count confound is real and quantified:** subject 1 runs 124 bpm at low load and 188
bpm at high load — **1.52× more beats per unit time**. A fixed-duration window therefore contains
far more samples at high power, and most irreversibility estimators are sample-size biased.
**Match sample counts, not durations.** This alone could manufacture the entire effect.

### Design (full version in `TODO.md` §C1)

1. Segment by recorded power into stages; take stationary windows within each.
2. Time irreversibility as EPR proxy, **≥2 independent estimators** — increment-asymmetry
   indices, plus a compression-based forward-vs-reversed KL estimator (Roldán–Parrondo spirit,
   which also connects to §5).
3. **Surrogate nulls:** shuffled and phase-randomised series must return zero. Without this the
   measure is uninterpretable.
4. Match sample counts across stages (see confound above).
5. Normalise power by each subject's ventilatory threshold so subjects are comparable.
6. Fit I ∝ P^n per subject, then pool.

**Read-out:** n ≈ 2 → tier II, same class as the colloids on this axis. n ≠ 2 → beyond linear
response, the first positive evidence that living systems need a different principle. Compare
against our measured colloid exponent (2.00, flat to 4.8%).

**Limitations to state up front:** mechanical power is a *compound* drive (metabolic throughput
*and* autonomic asymmetry), not a clean thermodynamic force — the *shape* of the scaling is the
observable, not the coefficient. RR irreversibility is a proxy for EPR, not EPR. Non-stationarity
within stages. 18 adolescent athletes is a narrow population. **A negative result (n ≈ 2) is
genuinely informative and should be reported as such.**
