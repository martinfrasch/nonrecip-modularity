# Where this system sits in the family of variational principles — and what that implies for living systems

Note written 2026-08-30, connecting the colloid results to the NWAP programme
(`minAction.net`) and the vertical/horizontal-scales manuscript.

---

## 1. Onsager holds here, and holds hard

Entropy production per particle, paired-bias-subtracted, at fixed configuration:

| χ | 1.5 | 2 | 3 | 5 | 8 |
|---|---:|---:|---:|---:|---:|
| EPR | 6.53e-3 | 1.22e-2 | 2.63e-2 | 7.25e-2 | 1.88e-1 |
| **EPR/χ²** | **2.90e-3** | **3.04e-3** | **2.93e-3** | **2.90e-3** | **2.93e-3** |

**EPR = kχ² with k ≈ 2.93×10⁻³, flat to 4.8% across a fourfold range of χ and a fifteenfold
range of EPR.** This is a membership test that could have failed at four independent points and
did not. The system is in the linear-response (Onsager) tier across its entire accessible range,
and a Rayleighian — minimise Φ̇ + Ψ with Ψ quadratic in the rates — is the appropriate
variational object for its dynamics.

## 2. The more interesting result is negative

This system has **nonreciprocity, sustained dissipation, spontaneous structure formation,
long-lived organised states, and a nucleation-separated two-phase architecture.** Every
ingredient usually invoked as a hallmark of active/living matter. And it is *still* in the
Onsager tier.

**So nonreciprocity + dissipation + self-organisation is not sufficient to leave linear
response.** Whatever distinguishes living organisation, it is not merely being driven and
non-reciprocal. That is a sharp constraint, and it makes this system valuable as a **control
for "driven but not alive"** — a well-characterised negative case against which candidate
biological signatures can be tested.

## 3. The family, and where the Triple-Action sits

| tier | regime | variational object | signature |
|---|---|---|---|
| I | equilibrium | free energy F = E − TS; MaxEnt | EPR = 0; detailed balance; Boltzmann |
| II | linear response | Rayleighian R = Φ̇ + Ψ, Ψ quadratic | **EPR ∝ drive²; Onsager reciprocity; single T_eff** |
| III | far from equilibrium | no general principle; **Maximum Caliber** is the leading candidate | EPR non-quadratic; reciprocity fails; observable-dependent T_eff |

Jaynes' **Maximum Caliber** — maximise path entropy subject to dynamical constraints — is the
natural home for a framework like NWAP, for a structural reason: it *contains an information
term by construction*, and it reduces to MaxEnt statically and to Onsager near equilibrium.
NWAP's Triple-Action (I_max − E_min + symmetry) has exactly that shape: an information term
traded against an energy term under structural constraints.

That is a placement worth pursuing, because it converts NWAP from a bespoke functional into a
**constrained member of a principled family** — one whose constraints are network-structural
rather than thermodynamic. It also predicts the limits: NWAP should reduce to a Rayleighian
when the drive is weak, and to free-energy minimisation when it vanishes. **Both reductions are
checkable**, and failure of either would be informative.

## 4. What keeps this falsifiable

The risk in "different systems need different principles" is that it explains everything and
forbids nothing. The discipline is that **each tier carries its own experimental signature, and
membership must be measured**:

1. **Quadratic dissipation** — EPR ∝ drive² at fixed configuration. *Run here; passed decisively.*
2. **Onsager reciprocity** — L_ij = L_ji among cross-responses. *Not yet run; requires a second
   independent drive. The sharpest remaining test.*
3. **Single effective temperature** — one T_eff reconciles response with fluctuation across
   observables. *Not yet run.*

A system passing all three is tier II and the Rayleighian applies. A system failing them is not,
and claiming a variational principle for it requires new evidence, not reassignment.

**The honest boundary from our own data:** EPR is quadratic *at fixed structure*, but the
structural response to χ is strongly non-linear — cluster count is non-monotonic, dipping at
χ=0.25 before rising tenfold. So the Rayleighian governs the dynamics at fixed structure;
**structure selection is a non-equilibrium transition that no current variational principle
covers.** Naming that gap is more useful than papering over it.

---

## 5. Vertical principles and horizontal scales: a concrete method from this work

The other manuscript asks how to identify horizontal scales of organisation from a vertical
(cross-scale) ansatz. This project produced a method, by accident, and a caution.

### The method: classify scale fixed points by stability

We asked what cluster size the system selects, and measured it operationally — per-cluster
**split rate** and **merge rate** as functions of size S, resolved at a timestep short enough
to see individual events:

| S | ~2 | ~5 | ~10 | ~19 | ~45 | ~74 |
|---|---:|---:|---:|---:|---:|---:|
| split rate | 0.333 | **0.453** | 0.254 | 0.181 | 0.081 | 0.086 |
| merge rate | 0.319 | 0.312 | **0.334** | 0.283 | 0.269 | 0.145 |

The rates cross at S\* ≈ 7–8. **Crucially, the crossing is a repeller, not an attractor:**
splitting dominates below, merging dominates above, so clusters flee S\* in both directions.

That distinction is the transferable content:

> **An attractor in the scale-flow is a populated level of organisation. A repeller is a
> boundary between levels.**

A system's horizontal scale ladder is therefore read off from the *fixed points of the
aggregation/fragmentation flow, classified by stability*. Attractors tell you where matter
accumulates; repellers tell you where one level ends and the next begins. Our system has a
repeller at S\*≈8, which creates two levels — a subcritical vapour below and an unbounded
condensate above — while populating neither at S\* itself.

### The operational recipe

1. Identify the scale variable S (cluster size, module size, community size).
2. Measure the *rates* of upward flow (merge/growth) and downward flow (split/fragmentation) as
   functions of S. **Rates, not net drift** — see the caution below.
3. Find crossings; classify each by the sign of (merge − split) on either side.
4. Attractors → levels. Repellers → level boundaries. No crossing → no scale selection at all.

### The caution, learned the hard way

Our first attempt used net drift ⟨dS⟩ over the production snapshot interval and returned +280
to +1284 in *every* size bin — nonsense, because over that interval small clusters are absorbed
into the condensate and the measurement reports condensate size rather than a growth increment.
The same error class recurred three times in this project under different names (σ²_v as a
"velocity variance", edge turnover as a "current", cluster displacement as a "speed").

> **A coarse-grained increment is not a rate.** Any cross-scale method that reads scale
> structure off differences taken at the wrong interval will manufacture attractors and
> boundaries that do not exist.

This is worth stating explicitly in a methods paper on identifying horizontal scales, because
the failure is silent — it produces plausible numbers with the wrong sign structure.

## 6. A testable programme for the living-systems question

The tier framework converts "which principle governs living organisation?" into a measurement:

- **Prediction:** biological systems should show **non-quadratic EPR scaling** in their drive,
  placing them outside tier II. If a driven biological system turned out quadratic, it would be
  in the same class as this colloid suspension — organised, dissipative, and unremarkable.
- **Control:** this colloid system is the negative case. Nonreciprocal, dissipative,
  self-organising, and still Onsager. Any claimed biological signature must distinguish itself
  from *this*, not merely from equilibrium.
- **Next measurement, either domain:** Onsager reciprocity (test 2 above). It is sharp, standard,
  and to our knowledge unrun on a nonreciprocal colloidal system. Passing places a system firmly
  in tier II; failing is the first positive evidence for needing something beyond it.

The exciting version of the question is therefore not "does a variational principle exist for
living systems" but **"at which tier do they sit, and does anything biological fail the tier-II
tests that this colloid passes?"** That is answerable.
