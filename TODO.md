# Systematic programme

Standing objectives, ordered. Status as of 2026-08-30.

---

## A. Tier-membership tests (which variational principle governs this system)

The discipline: a framework that says "some variational principle applies" forbids nothing.
Each tier carries a distinct experimental signature, and membership must be **measured**.

| # | test | signature of tier II (Onsager) | status |
|---|---|---|---|
| **A1** | **Quadratic dissipation** | EPR ∝ (drive)² at fixed configuration | ✅ **PASSED** — EPR = kχ², k ≈ 2.93e-3, flat to **4.8%** across χ = 1.5–8 (4× drive range, 15× EPR range) |
| **A2** | **Onsager reciprocity** | L_ij = L_ji among cross-responses | ⬜ **not run — sharpest remaining test** |
| **A3** | **Single effective temperature** | one T_eff reconciles response with fluctuation across observables | ⬜ not run |

### A2 — Onsager reciprocity (highest priority)

Requires a **second independent drive**; χ alone only varies one. Neither T1 (χ magnitude) nor
T3 (density, a static parameter) provides one.

- **Implementation:** uniform external force *f* applied to type-I particles only. Physically a
  field or gravity acting on one species — also a real experimental knob.
- **Fluxes:** J₁ = drift velocity of type-I relative to type-II (conjugate to *f*);
  J₂ = entropy production or cluster current (conjugate to χ).
- **Test:** measure the 2×2 response matrix; check ∂J₁/∂χ = ∂J₂/∂f after normalising by the
  thermodynamic forces.
- **The hard part is not the simulation.** It is identifying the conjugate force–flux pairs
  correctly; get that wrong and the symmetry test is vacuous. Think this through before writing
  code.
- To our knowledge this has not been run on any nonreciprocal colloidal system.

### A3 — effective temperature

Measure the response to a small perturbation and the corresponding fluctuation; test whether a
single T_eff reconciles them across observables. One T_eff → tier II. Observable-dependent
T_eff → tier III.

### A4 — the boundary our own data already mark

EPR is quadratic **at fixed structure**, but the structural response to χ is strongly
**non-linear** (cluster count non-monotonic: dips at χ=0.25, then rises tenfold).

> The Rayleighian governs the dynamics at fixed structure. **Structure selection is a
> non-equilibrium transition that no current variational principle covers.**

That gap is where new theory would have to go. Naming it is more useful than papering over it.

---

## B. Horizontal scale structure selection

**Goal:** predict the optimal organisation at each horizontal scale, and which elements comprise
each scale, from a vertical (cross-scale) principle.

### B1 — the method this project produced

Measure the **rates** of upward flow (merge/growth) and downward flow (split/fragmentation) as
functions of the scale variable S, then classify the crossings:

> **An attractor in the scale-flow is a populated level of organisation.
> A repeller is a boundary between levels.**

Demonstrated here: split and merge rates cross at S\* ≈ 7–8 particles, and the crossing is a
**repeller** — creating two levels (subcritical vapour, unbounded condensate) while populating
neither at S\* itself. Recipe in `VARIATIONAL_FAMILY.md` §5.

**Methodological caution, learned three times the hard way:** a coarse-grained increment is not
a rate. Net drift over the production interval returned +280 to +1284 in *every* size bin —
nonsense produced silently, with plausible magnitudes and the wrong sign structure.

### B2 — the renormalisation-group connection (to develop)

The scale-to-scale map with fixed points classified by stability **is** an RG flow. This gives
the "decompression" intuition a rigorous home: unfolding structure from a compressed
specification is an inverse RG flow, and the horizontal scales are its fixed points.

- **RG attractor (IR)** = stable scale, organisation accumulates → a populated level
- **RG repeller (UV / critical)** = boundary between phases → a level separator
- Relevant/irrelevant operator classification ↔ which perturbations survive coarse-graining ↔
  which elements comprise a given scale

**To do:** formalise the split/merge rate flow as an explicit RG recursion and check whether its
fixed-point structure reproduces the measured S\*. If it does, the method generalises beyond
cluster size to any scale variable.

### B3 — extracting the falsifiable core of the decompression picture

The metaphysical framing (an initial capacity computing optimal solutions from subquark to
cosmological scales) is not itself testable. The **structural claim inside it is**:

> If organisation at scale n+1 is decompressed from scale n, then scale n+1 must be largely
> predictable from scale n by a scale-local rule.

Measurable as **conditional description length / predictive information across adjacent scales**:

- H(scale n+1 | scale n) ≪ H(scale n+1) → decompression-like generation
- H(scale n+1 | scale n) ≈ H(scale n+1) → scales are independent; no decompression
- Estimable by compression-based methods (Lempel-Ziv, context-tree weighting) on coarse-grained
  representations at successive scales
- **Null model required:** the same measurement on a system with matched marginal statistics but
  scrambled cross-scale structure. Without that null the quantity always looks impressive.

This converts a vision into a number with a null. Keep the vision; publish the number.

### B4 — the vagus as a vertically organising substrate

The postulate: the vagus couples organisational levels (cellular → organ → organism → behaviour).

**Testable form:** if the vagus is a vertical organiser, cross-scale directed information should
be concentrated through it — so transfer entropy between scale levels should drop sharply under
vagal blockade or vagotomy, while within-level dynamics persist.

- Requires paired intact/blocked recordings with simultaneous multi-scale observables.
- The measurement is directed and conditional (TE between coarse partitions, conditioned on
  candidate blanket variables), not correlational.
- Connects to A2/A3: if the vagus carries the cross-scale flux, it is a natural place to look
  for a second drive and hence for reciprocity testing in a biological system.

---

## C. Biological vs non-biological: a concrete test with open data

**The question, made answerable:** not "does a variational principle exist for living systems"
but **"at which tier do they sit, and does anything biological fail the tier-II tests this
colloid passes?"**

**Why this is sharp:** our colloid suspension is nonreciprocal, dissipative, self-organising,
with long-lived structured states — and it is *still* Onsager (A1). So it is a **control for
"driven but not alive."** Any claimed biological signature must distinguish itself from this,
not merely from equilibrium.

### C1 — the experiment

**Prediction:** biological systems show **non-quadratic** scaling of entropy production with
drive, placing them outside tier II. If a driven biological system came out quadratic, it would
be in the same class as this colloid suspension — organised, dissipative, and unremarkable.

**Dataset:** [Cardiorespiratory measurement from graded cycloergometer exercise testing
(ACTES), PhysioNet v1.0.0](https://physionet.org/content/actes-cycloergometer-exercise/1.0.0/)
— 18 maximal graded exercise tests, beat-to-beat RR intervals, oxygen consumption **and
mechanical power output**, with ventilatory thresholds 1 and 2 per subject.

This is a close structural match to what we did: a **graded, directly measured drive** with a
simultaneous beat-to-beat response, exactly parallel to EPR(χ).

**Design:**

1. Segment each test into workload stages using the recorded power output.
2. Within each stage, take stationary windows of RR intervals.
3. Compute **time irreversibility** as an EPR proxy, with **at least two independent
   estimators** (increment-asymmetry indices, and a compression-based forward-vs-reversed KL
   estimator in the Roldán–Parrondo spirit — which also connects directly to B3).
4. Normalise against **surrogate nulls**: time-shuffled and phase-randomised series must give
   zero irreversibility. Without this the measure is uninterpretable.
5. Control for the confound that mean RR falls with load, so a fixed-duration window contains
   more beats at high power — match sample counts, not durations.
6. Fit irreversibility I ∝ P^n across stages, per subject, then pool.
7. Normalise power by each subject's ventilatory threshold so subjects are comparable.

**Read-out:**
- n ≈ 2 → tier II; physiology in the same class as the colloids on this axis
- n ≠ 2 → beyond linear response; the first positive evidence that living systems need a
  different variational principle
- Compare directly against our measured colloid exponent (2.00, flat to 4.8%)

**Stated limitations, up front:**
- Mechanical power is a **compound** drive (metabolic throughput *and* autonomic asymmetry), not
  a clean thermodynamic force. The *shape* of the scaling is the observable; the absolute
  coefficient is not comparable to the colloid k.
- RR irreversibility is a **proxy** for entropy production, not entropy production.
- Non-stationarity within stages; 18 subjects; adolescent athletes — a narrow population.
- A negative result (n ≈ 2) would be genuinely informative and should be reported as such.

### C2 — follow-on if C1 shows n ≠ 2

Then run the biological analogues of A2 and A3: a second drive (posture, temperature, or
pharmacological autonomic block) for reciprocity, and response-vs-fluctuation for effective
temperature. B4's vagal blockade design supplies a natural second drive.

---

## D. Housekeeping

- ⬜ Split the network-observable critique out of `PAPER.md` into a standalone methods paper
  (modularity degeneracy, turnover contamination, absent circulation) — see `FOLLOWUP.md` §3.
- ⬜ Add the coarsening law z(χ) and t_x ~ N^(1/z) as the closing result of the main manuscript.
- 🔄 T3 density continuations running (L=38, L=48 to t̂=1.6×10⁶) — settle whether the dilute
  many-cluster state is genuinely stable or merely below the binodal.
- ⬜ E1/A2 second-drive implementation in `simulate.py`.
