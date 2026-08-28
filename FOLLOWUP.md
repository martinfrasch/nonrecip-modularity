# Follow-up programme: timescales for hypothesis 1, and the variational family

Written 2026-08-27, after the transience result. Part 1 answers "on what timescale can
hypothesis 1 be affirmed?" with a measurement already in hand plus three designed experiments.
Part 2 addresses whether NWAP must broaden from Lagrangians to a family of variational
principles, and how to do that without making it unfalsifiable. Part 3 splits the material
between this manuscript and successors.

---

# Part 1 — The timescale question, made quantitative

## What we can already say

The finite-cluster state is not a steady state, but it decays by a **slow power law**, and that
law is measurable from data already collected. Fitting the largest cluster as S_max(t) = A t^z
over the growth regime (0.05N < S_max < 0.5N), at N=4000:

| χ | z | A | 1/z | t_x(N=4×10³) | t_x(N=2.2×10⁴) | t_x(N=10⁶) |
|---|---:|---:|---:|---:|---:|---:|
| 0 | **0.089** | 390 | 11.25 | 9.6e7 | 2.1e16 | 9.1e34 |
| 0.25 | 0.418 | 25.8 | 2.39 | 3.3e4 | 1.9e6 | 1.8e10 |
| 0.5 | 0.405 | 23.4 | 2.47 | 5.9e4 | 4.0e6 | 5.0e10 |
| 0.75 | 0.367 | 28.4 | 2.73 | 1.1e5 | 1.1e7 | 3.8e11 |
| 1.0 | 0.272 | 36.7 | 3.68 | 2.5e6 | 1.3e9 | 1.6e15 |
| 1.5 | 0.371 | 26.1 | 2.69 | 1.2e5 | 1.2e7 | 3.4e11 |

t_x is the crossover time at which the largest cluster reaches half the system — the end of the
finite-cluster regime.

**The answer to the timescale question: hypothesis 1 holds for t < t_x, and t_x ~ N^(1/z) with
1/z between 2.4 and 3.7.** A tenfold larger system stays in the finite-cluster regime 250 to
5000 times longer. In a macroscopic suspension the regime is effectively permanent, and the
distinction between "transient" and "steady state" becomes operationally empty.

That is a defensible reconciliation: the source paper is right about its system, we are right
that the state is not asymptotic, and the crossover scales as roughly the cube of system size.

## An inversion worth noting

χ=0 has by far the *slowest* coarsening (z = 0.089, near-arrest), and every nonreciprocal case
sits near the diffusion-limited value 1/3. **Nonreciprocity accelerates growth of the majority
phase** — it unjams the kinetic gel (`COARSENING.md`) — while simultaneously sustaining a
fragment population.

So "arrested coarsening" in this model means *a persistent population of small clusters*, not a
slowed majority phase. The genuinely arrested state is the reciprocal one. This should be stated
explicitly, because the phrase invites the opposite reading.

## Experiment T1 — does genuine arrest exist at stronger nonreciprocity?

We only ever tested χ ≤ 1.5. If z(χ) → 0 above some χ*, arrest becomes real rather than slow,
and hypothesis 1 is affirmed asymptotically in that regime.

- Run χ ∈ {2, 3, 5, 8} at N=4000, t̂=8×10⁵, 3 seeds. Cost ≈ 12 jobs × ~11 h.
- Measure z(χ) and lcf at convergence.
- **Prediction to be falsified:** z decreases monotonically above χ≈1 and extrapolates to zero at
  finite χ*. If instead z stays near 1/3 at all χ, no arrest transition exists and hypothesis 1
  is a finite-time statement at every coupling strength.
- χ > 1 is outside the derived EHD model, so this is a statement about the minimal model, not
  the colloid system. Say so.

## Experiment T2 — RESULT: the characteristic size is a critical nucleus, not a stable scale

**Run and answered.** The first attempt failed on existing trajectories: at the production
snapshot interval (Δt=2000) a small cluster is usually absorbed into the condensate within one
interval, so its plurality successor *is* the condensate and ⟨dS⟩ merely reports condensate size
(+280 to +1284 in every bin — nonsense). Repeating from equilibrated configurations with dense
sampling (Δt=50, 3 seeds, 105,708 cluster-observations) resolves individual merge and split
events:

| cluster size | n | split rate | merge rate | which dominates |
|---:|---:|---:|---:|---|
| ~2 | 63,978 | 0.3326 | 0.3190 | split (marginal) |
| ~5 | 24,978 | **0.4525** | 0.3118 | **split — clusters dissolve** |
| ~10 | 10,076 | 0.2542 | **0.3335** | **merge — clusters grow** |
| ~19 | 3,592 | 0.1812 | 0.2831 | merge |
| ~45 | 1,832 | 0.0813 | 0.2686 | merge |
| ~74 | 502 | 0.0857 | 0.1454 | merge |

**The rates cross at S\* ≈ 7–8 particles, and the crossing is a repeller, not an attractor.**
Below S\* splitting dominates and clusters dissolve; above S\* merging dominates and they grow
without bound. That is a **critical nucleus**, the classic nucleation signature — an *unstable*
fixed point.

This is the decisive answer for hypothesis 1. The E-vs-C tradeoff requires a **stable** attractor
at intermediate size, so that clusters larger than S\* are driven back down by connection cost.
The system has the opposite: anything above S\* grows without limit. There is no mechanism here
that caps cluster size.

It also explains T4 exactly. The fragment population (mean ≈ 5) sits **below** the critical
nucleus — these are subcritical clusters, continually forming and dissolving. Their size
statistics are N-independent because S\* is set by local energetics, not by system size. The
steady state is a condensate coexisting with a **subcritical vapour**, and the "characteristic
scale" is the nucleation barrier separating them.

Caveat: ⟨dS⟩ remains positive in every bin even at Δt=50, because plurality-tracking still
catches absorption events. The split/merge *rates* are the reliable measure; ⟨dS⟩ is not.

## Experiment T2 (original design, superseded by the result above)

A characteristic size S* exists if and only if the fission rate and coalescence rate cross at
finite size. This is the measurable form of the energy-versus-connection-cost tradeoff, and it
predicts S* without running to steady state.

- Track clusters between consecutive snapshots; count merge and split events resolved by cluster
  size. Existing trajectories suffice — no new simulation.
- Extract k_coal(S) and k_fiss(S); S* is where they cross.
- **Prediction:** at χ ≤ 1.5 they do not cross (coalescence dominates at all S, hence phase
  separation). Extrapolating k_fiss(S,χ)/k_coal(S,χ) in χ predicts the χ* at which they first
  cross — which **T1 then tests directly**. Two independent routes to the same number is the
  strongest available check.

## Experiment T3 — density, the other knob

Coalescence requires clusters to find each other, so dilution should extend the regime
independently of χ. At fixed N, vary the box to span packing fractions 10–40%.

- **Prediction:** t_x grows steeply as density falls, and below some density finite clusters are
  stable for any accessible time. If a genuine dilute-limit arrest exists, that is where
  hypothesis 1 is affirmed in the strongest sense.
- This is also the most experimentally actionable prediction, since colloid density is the
  easiest parameter to vary in the lab.

## Experiment T4 — RESULT: yes, and it is the subcritical vapour

**Run and answered.** Fragment statistics excluding the largest cluster, late window:

| N | mean fragment | median | p95 | p99 | fragment mass fraction |
|---:|---:|---:|---:|---:|---:|
| 4,000 | 5.02 | 3 | 14 | 33 | 17.6% |
| 9,000 | 4.78 | 3 | 12 | 35 | 15.8% |
| 16,000 | 5.34 | 3 | 11 | 26 | 16.9% |
| 22,000 (source spec, χ=1) | 7.15 | 3 | 19 | 53 | 14.7% |

Mean fragment size is constant at ≈5 particles across a **5.5× range of system size**, with an
N-independent cutoff and a constant ~16% mass fraction. **A scale is selected** — the earlier
statement that no scale exists anywhere was too strong.

But the distribution is monotonically decaying (11474, 4109, 1225, 314, 76, 28, 25 by octave),
not peaked, and T2 identifies the scale as the critical nucleus. So the selected scale is a
nucleation barrier, not a modular cluster size. Hypothesis 1 survives only in the weak sense
that *something* is scale-selected; the mesoscale attractor its tradeoff requires does not exist.

A secondary test — whether fragment mass fraction rises with χ, as a saturation vapour pressure
would — is **inconclusive**: 0.212, 0.096, 0.195, 0.457, 0.203 for χ = 0.25…1.5, with errors up
to ±0.288. Fluctuations dominate at 3 seeds.

## Experiment T4 (original design, superseded by the result above)

We showed the *condensate* has no characteristic size. We have not asked whether the fragments
do. Median fragment size is 3 at every box and every χ, which is suspiciously stable.

- Measure the fragment size distribution excluding the largest cluster, at all box sizes.
- **Prediction:** the fragment distribution is scale-free (power law) with an exponential cutoff
  set by the fission mechanism. If instead it is peaked at a finite size independent of N, then
  **a characteristic scale does exist — in the fragments, not the condensate** — and hypothesis 1
  would be affirmed in a modified form that nothing so far has tested.

This is the cheapest experiment and the one most likely to rescue hypothesis 1.

---

# Part 2 — Does NWAP need a family of variational principles?

## Yes, but the family already exists and is not ad hoc

The important point is that this is **not** a broadening invented to accommodate an
inconvenient result. There is a standard hierarchy of variational principles for stochastic
dynamics, and our measurements place this system in a specific tier of it:

| tier | condition | variational object | signature |
|---|---|---|---|
| **I. Equilibrium** | F = −∇U | free energy; stationary state minimises it | EPR = 0 exactly; detailed balance; Boltzmann |
| **II. Linear response / Onsager** | F = −∇U + v, v small | Rayleighian R = Φ̇ + Ψ, Ψ quadratic in rates | **EPR ∝ (drive)²**; Onsager reciprocity L_ij = L_ji; FDT with an effective temperature |
| **III. Far from equilibrium** | v large | no general principle; Onsager–Machlup action still defined but not a selection principle | EPR non-quadratic; reciprocity fails |

Crucially, **the Onsager–Machlup action exists for any drift field**, including non-gradient
ones: S[x] = ∫ (ẋ − μF)²/4D dt. Least action does not fail for nonreciprocal systems. What fails
is the *equilibrium* corollary — that the stationary distribution is exp(−βU) and structure is
selected by minimising an energy. That corollary is what the original NWAP framing relied on.

## Our data place this system in tier II, and that placement was measured

At fixed configuration, EPR = kχ² (net/χ² flat to 1.27× and 1.61× on two of three
configurations, `FINAL_RESULTS.md` §4). That is the tier-II signature, measured rather than
assumed. Had it come out cubic, the Rayleighian would not apply and we would have had to say so.

## How to broaden without becoming unfalsifiable

The user's worry is the right one. The discipline that prevents it:

**Every tier carries its own experimental signature, and membership must be measured, not
asserted.** Concretely, three tests any candidate system must pass to claim tier II:

1. **Quadratic dissipation.** EPR ∝ (drive)² at fixed configuration. *Falsifiable, and we ran it.*
2. **Onsager reciprocity.** Apply two independent perturbations, measure the cross-response
   coefficients, check L_ij = L_ji. *Not yet run here — see E1 below. This is the sharpest
   available test and the one that would most strengthen the framework.*
3. **Fluctuation–dissipation with a single effective temperature.** Compare the measured
   response to a small perturbation against the equilibrium-form prediction from the
   fluctuations. *Not yet run — see E2.*

A framework that says "some variational principle applies" is unfalsifiable. A framework that
says "this system is tier II, therefore EPR is quadratic, reciprocity holds, and FDT holds with
one effective temperature" makes three independent predictions that can each fail. That is the
version worth writing.

## The honest boundary of the claim

There is an important limit visible in our own data. EPR is quadratic **at fixed configuration**,
but the structural response to χ is strongly **non-linear** — cluster count is non-monotonic,
dipping at χ=0.25 before rising tenfold. So:

> The Rayleighian applies to the dynamics *at fixed structure*. Structure selection is a
> non-equilibrium transition that no current variational principle covers.

Stating that boundary explicitly is what separates a framework from an unfalsifiable one. It also
identifies exactly where new theory would be needed, which is more useful than claiming coverage.

## Two experiments that would test the tier assignment

**E1 — Onsager reciprocity.** Introduce a second small drive (e.g. a weak uniform field, or a
second reciprocity parameter acting on a different pair type) and measure the 2×2 matrix of
cross-responses between the two drives and two conjugate fluxes. Check symmetry.
*This is the single most valuable next measurement for the theory paper*: it is a sharp,
standard, falsifiable test of tier-II membership that nobody has run on a nonreciprocal
colloidal system.

**E2 — Effective-temperature consistency.** Measure the response to a small perturbation and the
corresponding fluctuation, and test whether one effective temperature reconciles them across
observables. A single T_eff supports tier II; observable-dependent T_eff indicates tier III.

---

# Part 3 — What belongs in this manuscript, and what does not

## Keep in the current manuscript

It is already a coherent single story: **a method to isolate the antisymmetric sector, and what
that sector does in one well-characterised system.**

- the χ decomposition, its validation, and the demonstration that α̂ and the steric size ratio
  cannot serve as reciprocity knobs
- nonreciprocity as the causal driver of the cluster statistics (12.4×, p=9.4e-10)
- the reciprocal reference as a kinetic gel; unjamming; the non-monotonic n_cl
- entropy production, positive and quadratic
- box scaling, the transience result, and the source-specification replicate
- **add: the coarsening law z(χ) and the crossover time t_x ~ N^(1/z)** — this converts the
  transience result from a negative finding into a quantitative one, and it is the natural
  closing result of the paper
- **add: the inversion** — that χ=0 is the arrested case and nonreciprocity accelerates the
  majority phase. It corrects a natural misreading of "arrested coarsening" and costs one
  paragraph.

## Move out of this manuscript

- **The network-observable critique** (modularity degeneracy, turnover contamination, absent
  circulation) is currently a major theme and it dilutes the paper. It is a genuinely useful
  methods contribution for the network-physics community and stands better alone: *"Network
  observables do not detect nonreciprocity: modularity, partition turnover and circulation on
  contact graphs."* Keep a compressed version in the present paper as a methods caution with a
  forward reference.
- **The Markov-blanket material** is speculative and belongs nowhere near a physics paper on
  colloids.

## Next paper 1 — theory

*"Nonreciprocal networks and the antisymmetric sector: which variational principle survives."*
Part 2 above is its skeleton. Empirical anchor: EPR ∝ χ² from this work, plus E1 (reciprocity)
and E2 (effective temperature). Its contribution is the tier framework with falsifiable
membership tests, and an explicit statement of where variational coverage stops. Cite Hara et
al. and Fruchart et al. as physics context.

## Next paper 2 — kinetics

*"Arrest without a steady state: coarsening kinetics of a nonreciprocal colloidal suspension."*
Built on T1–T4. Its contribution is t_x ~ N^(1/z), the fission/coalescence balance, and either
the location of a genuine arrest transition at χ* or a demonstration that none exists. This is
where hypothesis 1 gets its fair and complete test.

## Recommended order

T4 and T2 first — both use existing trajectories, cost nothing, and either could change the
verdict on hypothesis 1 before any new compute is spent. Then T1 (arrest transition) and E1
(reciprocity), which are the two measurements that most affect the theory paper. T3 last: it is
the most experimentally relevant but the least likely to overturn anything.
