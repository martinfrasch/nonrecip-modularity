# Isolating the antisymmetric sector of a nonreciprocal colloidal model: kinetic unjamming, transient arrested coarsening, and irreversibility without coarse-grained circulation

**Martin G. Frasch**

*Computational study, August 2026. Repository: `nonrecip-modularity`.*

---

## Abstract

Nonreciprocal interactions do not derive from a scalar potential, which removes the equilibrium
apparatus — detailed balance, and with it structure selected by energy minimisation — even
though a path-space variational representation survives; a Boltzmann-shaped stationary density is
not excluded in general, only detailed balance is, and whether the density survives is a question
for each model. A widely held intuition holds that
such a coupling can be split into a symmetric sector that selects structure and an antisymmetric
sector that merely adds circulation without disturbing it. We test that intuition in the
agent-based model of Hara et al. [Phys. Rev. Lett. **137**, 068302 (2026)] for size-asymmetric
colloids driven by electrohydrodynamic flows.

We introduce a reciprocity mixing parameter χ that scales the antisymmetric part of the pair
coupling while leaving the symmetric part bit-for-bit unchanged, with χ = 0 restoring Newton's
third law exactly and χ = 1 recovering the published force law. Two parameters previously used for
this purpose are shown to be unsuitable: the coupling strength cancels exactly from the
antisymmetric-to-symmetric ratio, and the steric size ratio, as parameterised in the source model,
varies packing at fixed electrohydrodynamic radii.

Across 113 blocked-design simulation runs — typically three shared-seed realisations per principal
condition, so that runs at different χ share initial conditions and noise streams — extended by 25
continuations, at five system sizes from N = 10³ to 2.2 × 10⁴,
we find that the antisymmetric sector alters *static* structure by an order of magnitude, refuting
the solenoidal intuition. It produces the finite-time fragmented morphology conventionally called
arrested coarsening while simultaneously *unjamming* the reciprocal gel and accelerating growth of
the majority phase, so that cluster count depends non-monotonically on χ and the genuinely arrested
state is the reciprocal one. That fragmented state is a long-lived transient rather than a steady
state: the largest cluster grows in proportion to system size (N^1.00 ± 0.07), the time to
condense scales linearly with N, and no stable characteristic cluster size is selected. The scale
that is selected, S* ≈ 8–10 particles, is a critical-size crossover: the committor of the
fragment-size dynamics, converged in the observation horizon, is one half at S ≈ 8, and the
monomer-exchange drift changes sign at S ≈ 10. Excess dissipation over the reciprocal control,
measured with a configurational estimator that requires no trajectory and agrees with the
trajectory estimator to 5–10 % where both resolve, is unresolved for χ ≤ 0.25 and reaches a
quadratic plateau by χ = 0.75; the configurational decomposition cancels 87–91 % of the unopposed
antisymmetric contribution above that crossover, through the steric contact force, and all of it
within error below, through the pair attraction. The same crossover governs structure: the
fragment population is approximately stationary over the measured windows, with each channel pair
— shedding and reabsorption, fission and fusion, nucleation and dissolution — approximately
balanced, and the non-monotonic cluster count is accompanied by distinct drive dependences of
condensate shedding and per-fragment reabsorption, which switch on at different χ. Because the
nonreciprocal forces do not sum to zero, the suspension also acquires a net internal force and a
centre-of-mass drift, near-ballistic over short windows but finite-size in amplitude. Yet no
circulation is detectable in any coarse observable we examined, at either the system-averaged or
the single-cluster level. That is a statement about those projections rather than about the
dynamics, since a vanishing projected current does not imply detailed balance even though the
converse holds. Several standard network measures, Newman modularity among them, prove unable to
detect nonreciprocity at all.

Three auxiliary strands are reported more briefly. Of three commonly invoked near-equilibrium
diagnostics, two fail for reasons concerning the instrument rather than the system; the third,
measured on single-particle coordinates, gives a single effective temperature equal to the bath
temperature over a drive-dependent short-lag window, shared by both species, with a shared and
growing violation at longer lags.
Comparisons with beating *Chlamydomonas* axonemes, with a synthetic ensemble given an imposed
turning rate, and with published social-dominance matrices bound how far the circulation result
generalises — in particular, the axoneme contrast cannot be read as biological, since the synthetic
ensemble circulates more strongly still.

**Keywords:** nonreciprocal interactions, active matter, arrested coarsening, entropy production,
kinetic arrest, coarse-graining

---
## 1. Introduction

Nonreciprocal interactions, in which the force that *i* exerts on *j* is not the negative of the
force *j* exerts on *i*, do not derive from a scalar potential. They are now recognised as a
generic route to phase behaviour with no equilibrium counterpart [2]. It is worth being precise about
what this does and does not preclude, because the two are often conflated. It does *not* preclude a
variational representation of the dynamics: the Onsager–Machlup action is well defined for any
drift field, gradient or otherwise, and path-space action principles apply to dissipative and
stochastic systems generally. What it generally removes is the *equilibrium corollary* — detailed balance, and with it
structure selected by minimising an energy. It need not exclude a Boltzmann-shaped stationary
density in general, although Section 3.2 shows that the added drift of this model does not
preserve the reciprocal Boltzmann density.

The questions that are actually at stake are therefore narrower and sharper than "does an action
principle survive". They are: does the system admit an equilibrium potential governing its
stationary distribution; does the antisymmetric part of the interaction preserve a reference
stationary density; and is a near-equilibrium Onsager–Rayleigh construction justified? A
variational representation of trajectories is not the same object as a physical extremum principle
that selects stationary structure, and this paper is concerned with the latter.

The question has become concrete rather than merely formal. Hara, Sumino and colleagues recently
reported a controlled experimental realisation: polystyrene colloids of two radii, confined
between indium tin oxide electrodes under an AC field, develop electrohydrodynamic flows whose
strength scales steeply with particle radius [1]. Size-asymmetric pairs therefore experience
imbalanced attraction and spontaneously self-propel. The resulting suspension exhibits *arrested
coarsening*: clusters continuously fragment and reorganise rather than growing without bound as
they do in the monodisperse case. Accompanying agent-based simulations identify nonreciprocal pair
propulsion as the minimal ingredient for this behaviour.

The intuition we test is a general one, and predates any particular framework: that a nonreciprocal
coupling may be decomposed into a symmetric part generating the gradient, energy-like component of
the dynamics, and an antisymmetric part generating a *solenoidal*, circulating component,
identified with the nonreciprocal propulsion that Hara et al. isolate as their minimal ingredient.
Variants of this decomposition appear across the nonreciprocal-matter literature [15–18]. The
immediate motivation for this study was one specific formulation: a *directed-graph extension* of
the Network-Weighted Action Principle [9], which assigns the antisymmetric sector of a directed
coupling the circulating, nonequilibrium role and predicts that it should leave a sustained
modularity excess in the resulting contact network [30]. We should be precise about the status of
that extension, since it is what we set out to test: the action principle itself is published
[9, 30, 33], but the directed-graph extension is not, and the specific assignment tested here was
formulated by the present author in the course of setting up this study. It is therefore a
hypothesis of ours, not a claim drawn from the literature, and we test it as such. We return to it
in Section 4.7; the predictions tested below are consequences of the decomposition itself and do
not depend on any particular framework endorsing it.

Three consequences follow from that assignment, and all three are testable. We note at the outset
that the assignment involves a step that is not automatic: antisymmetry under exchange of particle
labels, violation of Newton's third law, divergence-freedom of a vector field, and circulation of
the configuration-space probability current are four distinct properties, and a reciprocity
decomposition of the pair force is not the same object as a Helmholtz–Hodge decomposition of the
probability-current field. An antisymmetric pair-force component need not be divergence-free in the
full 2N-dimensional configuration space, and even a divergence-free drift need not preserve a given
density: the relevant condition is ∇·[ρ_ss(**X**) **v**(**X**)] = 0, not ∇·**v** = 0. Our results
below bear on whether the identification holds empirically in this system; they do not establish
that it was ever justified formally.

First, because a solenoidal field is divergence-free, the antisymmetric sector should not alter any
static structural observable. Second, its entire signature should therefore appear in probability
currents, giving "frozen modularity, circulating partition". Third, the cluster-size distribution
and reorganisation rate should be set by the ratio of antisymmetric to symmetric coupling
strength. The appeal of the picture is that it would preserve an action principle: the symmetric
part selects structure by energy minimisation, and the antisymmetric part adds circulation on top
without disturbing the outcome.

We test all three predictions. The first two fail. The third is borne out in substance — cluster
statistics are strongly controlled by the antisymmetric-to-symmetric ratio — but not in the
monotonic form the picture supposes, and we set out in Section 4.7 what that distinction costs it.

Testing that proposal requires a control parameter that varies the ratio of antisymmetric to
symmetric coupling and nothing else. We show in Section 3.1 that neither parameter previously
proposed for the role does so, and we construct one that does. With that instrument in hand, three
distinct questions become answerable, and they organise the paper.

The first is causal: what does nonreciprocity actually do to this system? Section 3.2 establishes
that it drives the fragmented morphology called arrested coarsening, and that it does so by
altering equal-time structure — which directly contradicts the solenoidal premise. Sections 3.3 to 3.7 then examine the observables
themselves, and find that several standard network measures are unable to detect nonreciprocity
at all. Sections 3.5 and 3.8 to 3.9 characterise the resulting steady state, showing that the
reciprocal reference is a kinetically arrested gel rather than an equilibrium structure, that no
finite characteristic cluster size is selected, and that the scale the system does select is a
nucleation barrier rather than a stable organisational level.

The second question is which variational tier the system occupies. Rather than asserting that some
principle applies, we treat tier membership as an experimental matter with measurable signatures,
and test three of them: quadratic dissipation, Onsager reciprocity, and a single effective
temperature (Section 3.10). Two cannot decide it. The quadratic dissipation is near-tautological
at fixed structure; what it does establish is a crossover in χ, below which the excess is
unresolved and the configurational decomposition cancels all of the unopposed contribution within
error, and above which it cancels 87–91 % and the residue is quadratic. An apparent Onsager antisymmetry does not survive a
control in which the two response channels are made structurally dissimilar. The third, measured
on single-particle coordinates with a force pattern conjugate to the single-particle fluctuation,
does decide something: a single effective temperature equal to the bath temperature, shared by
both species, holds over a short-lag window that shrinks as the drive grows, and a shared, growing
violation replaces it at longer lags. The tested single-particle fluctuation–response relation thus has a measured short-lag
agreement window rather than a tier, and we regard the route to that statement as a result about the difficulty of
applying the taxonomy as much as about the system. Chasing the third diagnostic on species
coordinates also turned up a mechanism: because nonreciprocal forces do not sum to zero over the system, there is a
net internal force and the centre of mass drifts, near-ballistically over short windows
(⟨ΔR²⟩ ~ t^1.96) where the reciprocal limit is diffusive (t^1.01). That accounts for about half of
the anomalous growth; the rest survives a control that removes it. The drift is finite-size —
its amplitude falls roughly as N^−1/2, as the incoherent addition of pair propulsions predicts —
and window-limited, decorrelating over long runs.

The third question follows from the second, though it does not depend on how the second is
resolved. Whatever tier the colloid suspension occupies, it is nonreciprocal, dissipative and
strongly self-organising while being uncontroversially not alive, which makes it a useful control:
any proposed signature of biological organisation must distinguish itself from this system, not
merely from equilibrium. Section 3.11 therefore uses the colloid suspension as such a control —
a well characterised system that is driven but not alive — and compares it against published
recordings of beating cilia. The two differ sharply, and on an axis that neither the modularity
measures nor the entropy production alone would have identified; a synthetic control then shows
that the axis is not the biological one it appears to be. Section 3.12 closes with a system in
which the symmetric/antisymmetric decomposition can be computed exactly rather than inferred —
published social-dominance matrices — and reaches the paper's central conclusion about antisymmetry
by that independent route.

We use the source model without modification except where explicitly stated, and we validate our
reimplementation against the published specification and duration in Section 3.8.

## 2. Model and methods

### 2.1 Equations of motion

We reimplemented the model of ref. [1] (End Matter Sec. II, Eqs. 7–9). Overdamped dynamics in
two dimensions with periodic boundaries, in nondimensional units (length λ = 9 μm, time
τ = 0.01 s):

    dx_i/dt = ξ_i(t) + (1/s_i) [ Σ_j F^col_ij + Σ_j F^EHD_ij ]

$F^{\rm col}$ is a reciprocal soft-core linear spring of natural length s_i + s_j. The EHD term is

    F^EHD_{i←j} = −α · l_j⁴ / (r² + l_j²)^{5/2} · **r**_ij ,    r ≤ 1

and is nonreciprocal because the force on *i* scales with *l_j*, the **other** particle's EHD
radius; the arrow in the subscript records that convention, which Section 2.2 makes explicit.
Noise satisfies ⟨ξ_i ξ_j⟩ = (σ²/s_i) δ_ij δ(t−t′), σ = 2.6×10⁻³. Equivalently, the mobility
tensor is μ_ij = s_i⁻¹ δ_ij **I** and the bare diffusion coefficient D_i = σ²/(2 s_i), so that
D_i/μ_i = σ²/2 for every particle irrespective of size. Baseline parameters
follow ref. [1]: α̂ = 0.005, l_L = s_L = 1/6, l_S = s_S = 1/9, steric polydispersity
ε ~ N(0, 1/30), composition N_L : N_S = 1 : 3.

Integration is Euler–Maruyama with dt̂ = 0.05 (stability bound ≈ 0.2), using a hand-rolled cell
list in a single compiled kernel.

### 2.2 The reciprocity mixing parameter χ

Write g(r; l) ≡ α l⁴/(r²+l²)^{5/2} for the electrohydrodynamic kernel, and define the
coefficients acting *on* each particle,

    a_i ≡ g(r; l_j) ,   a_j ≡ g(r; l_i)

**The index convention matters and is the whole content of the model:** the coefficient acting on
a particle is set by the *other* particle's EHD radius, so a_i carries l_j and a_j carries l_i. In
this notation the published law of ref. [1] reads F_i = −a_i **r**_ij, F_j = +a_j **r**_ij, and the
imbalance a_i ≠ a_j is the nonreciprocity. With ḡ ≡ (a_i+a_j)/2 we replace those coefficients by

    c_i = (1−χ)·ḡ + χ·a_i ,   c_j = (1−χ)·ḡ + χ·a_j
    F_i = −c_i **r**_ij ,     F_j = +c_j **r**_ij

so that the **symmetric part (c_i+c_j)/2 = ḡ is independent of χ**, while the antisymmetric
part (c_j−c_i)/2 = χ(a_j−a_i)/2 is exactly linear in χ. At χ=1, c_i = a_i and the law of ref. [1]
is recovered exactly; at χ=0, c_i = c_j = ḡ and Newton's third law holds identically. In the
kernel (`simulate.py`) these coefficients appear as `gi` and `gj`, with `gi` built from `l4[j]`;
the same convention is used in `fdt.py`, `onsager.py` and `simulate_rect.py`.

Validation against the compiled kernel (`tests/test_reciprocity.py`): an isolated noise-free
pair of equal steric radii, hence equal mobilities (so that zero net force implies a fixed
arithmetic centre; with unequal mobilities a reciprocal pair's arithmetic centre drifts as
(μ_i − μ_j)**F**_i/2 and only the friction-weighted centre is fixed), shows centre-of-mass drift
0.000×10⁰ at χ=0, 3.95×10⁻² at χ=1, and exactly half that at
χ=0.5 (1.975×10⁻², linear in χ to four digits); an equal-radius pair follows identical
trajectories at χ=0 and χ=1 to 10⁻¹⁴. Against the pre-χ code at χ=1, trajectories are bitwise
identical for 1000 steps, diverging only to 9.8×10⁻¹⁵ by 5000 steps (floating-point rounding
amplified by chaos).

### 2.3 χ=0 gives equilibrium-compatible dynamics

At χ=0 all forces are central, pairwise and equal-and-opposite, hence conservative; mobility
μ_i = 1/s_i with noise variance σ²/s_i gives D_i/μ_i = σ²/2, uniform across particles.
Fluctuation–dissipation therefore holds, detailed balance holds for the χ = 0 *equations of
motion*, and the stationary probability current vanishes.

A distinction must be maintained throughout, because we later show that it matters. The χ = 0
*dynamics* are equilibrium-compatible and possess a Boltzmann stationary distribution. The χ = 0
*configurations realised in finite-duration simulations* are kinetically trapped and demonstrably
do not sample that distribution (Section 3.5). We use "equilibrium-compatible reciprocal dynamics"
for the former and "kinetically trapped finite-time ensemble" for the latter, and reserve
"equilibrium" for the formal long-time distribution. The value of χ = 0 as a reference lies in the
first sense: it guarantees that any measured current at χ > 0 is attributable to the antisymmetric
coupling and not to the reference itself. This is a stronger reference than the monodisperse comparison of the
source experiment, which differs from the bidisperse system in reciprocity *and* polydispersity
simultaneously.

### 2.4 Observables

Contact graph: edge if r < 1.1(s_i+s_j), minimum image. Per snapshot we compute Newman
modularity Q [3] of the Louvain partition, modularity of a degree-preserving (double-edge-swap)
null, the number of connected components with ≥2 particles (n_cl), and the largest-cluster
fraction (lcf). Per consecutive pair we compute the adjusted Rand index between partitions
(reported alongside a same-graph two-seed noise floor), the Jaccard distance between contact
edge sets, and σ²_v.

Two further observables were introduced during this work:

**Circulation.** The signed area rate in a plane of two observables, A = ½⟨x ẏ − y ẋ⟩. At
stationarity, detailed balance implies A = 0 for *every* observable pair, so a nonzero value is
sufficient evidence that time-reversal symmetry is broken. **The converse does not hold**: A = 0
does not establish detailed balance, because an irreversible current may simply be absent from the
chosen two-dimensional projection. That asymmetry is not a technicality here — it is precisely what
Sections 3.4 and 3.11 report, and the reason the paper's negative circulation results are
statements about projections rather than about the dynamics.

**Dissipation.** The kernel accumulates the Stratonovich heat delivered to the medium,

    Q̇ = Σ_i **F**_i ∘ **ẋ**_i ,

by midpoint rule [5, 26]. The dynamics are overdamped with a uniform fluctuation–dissipation ratio
D_i/μ_i = σ²/2 ≡ T for every particle (Section 2.1), so the medium entropy-production rate is
Ṡ_med = Q̇/T. We report the per-particle rate Q̇/(T N); with entropy measured in units of k_B this
has dimensions of inverse time, per unit t̂ in the nondimensional units of Section 2.1. Under detailed balance ⟨Q̇⟩ = 0
identically, which is what makes this a current — unlike edge turnover or σ²_v, both of which are
nonzero at χ = 0.

Two caveats attach to the name. The trajectory estimator is unusable at the production timestep:
it subtracts two nearly cancelling terms of order μ|F|², and the discretisation residual dominates.
Where we use it, it is measured at dt = 6.25×10⁻⁴, restarting from equilibrated configurations,
and what it reports is a *difference*: each configuration is probed at its own χ and again at χ = 0
with the same noise stream, and the χ = 0 value subtracted. That removes the configuration-dependent
discretisation residual, which is driven by the symmetric forces and is common to both probes. This
is an effective numerical control, not a derivation, and it presumes the residual is χ-independent
at fixed configuration. We therefore speak of the **excess dissipation over the reciprocal
control**, and use "entropy production" only where the distinction does not bear on the argument.

**A configurational estimator of the excess.** The paired trajectory probe is not the only route.
Writing F = F_s + χF_a for the reciprocal and antisymmetric sectors and using the Stratonovich
conversion ⟨g∘ξ_i⟩ = D_i∇_i·g, the mean heat separates as ⟨Q̇⟩ = ⟨F_s∘ẋ⟩ + χ⟨F_a∘ẋ⟩ with
⟨F_s∘ẋ⟩ = −d⟨U⟩/dt, which is measured directly from the stored energies and is of order 10⁻⁶ per
particle in units of T per unit t̂ — the units of the dissipation tables below, whose entries are
of order 10⁻³ — in the late windows used. The remainder is a configurational average of known functions:

    χ ⟨F_a∘ẋ⟩ = χ Σ_i (1/s_i) ⟨F_a,i · F_i⟩ + χ Σ_i D_i ⟨∇_i · F_a,i⟩ ,

with F_a the pair-propulsion field, F the full force, and ∇·F_a evaluated analytically per pair
(plus a shell term for the cutoff at r = 1). The excess dissipation is therefore an average over
stored snapshots, with no trajectory, no small timestep and no subtraction of cancelling terms
(`static_epr.py`). It splits naturally as $\chi(\tilde b+\chi\tilde a)$, with $\tilde a$ = Σ⟨|F_a,i|²⟩/(s_i T N) ≥ 0 the
dissipation the antisymmetric force would produce unopposed and $\tilde b$ = [Σ⟨F_a·F_s⟩/s_i +
Σ D_i⟨∇·F_a⟩]/(T N) the part cancelled by the response of the contact forces; in a Boltzmann
ensemble $\tilde b$ vanishes identically by integration by parts, so $\tilde b$ ≠ 0 measures how far the structure
has rearranged under the drive. One correction is needed: the production timestep distorts
contact-pair statistics (dt·k·μ_max = 0.5), which shows up as a violation of the identity
⟨F_s∘ẋ⟩ = −dU/dt by 4.4 T per particle per unit time on stored snapshots. Re-equilibrating each
snapshot for t = 500 at dt = 0.005 reduces that residual to 0.35, and every configurational value
quoted below is taken from such re-equilibrated configurations (162 per run) unless stated. Where
both estimators resolve (χ = 1.5 to 8 at N = 4000), they agree to within 5–10 %.

### 2.5 Simulation campaign

**113 simulation runs, extended by 25 continuation runs, for 138 trajectory files in total**, at
five system sizes. Runs at different χ are *not* statistically independent: seeds are blocked, so
seed *k* gives the same initial configuration and noise stream at every χ. This is a
common-random-number design chosen to reduce the variance of χ contrasts, and it is why we quote
effect sizes rather than p-values in Section 3.2. Continuations extend existing runs and are not
replicates.

| series | N | box | t̂ | runs |
|---|---:|---|---:|---:|
| Pilot (baseline, α̂ sweep, size-ratio sweep, χ sweep) | 1,000 | L = 12 | 1×10⁵ | 61 |
| Paper scale (χ sweep 6 × 3, monodisperse reference × 3) | 4,000 | L = 24 | 4×10⁵ | 21 |
| Box scaling at fixed density 6.944 particles/area | 9,000 / 16,000 | L = 36 / 48 | 4×10⁵ | 6 |
| Density series at fixed N, χ = 1.5 | 4,000 | L = 30 / 38 / 48 | 8×10⁵ | 6 |
| High-drive series, χ ∈ {2, 3, 5, 8} | 4,000 | L = 24 | 8×10⁵ | 12 |
| Rectangular-box control, aspect 1.8 : 1 | 4,000 | 32.20 × 17.89 | 8×10⁵ | 3 |
| Replicate at the source specification, 22.7% type-I | 22,000 | L = 53.67 | 3.6×10⁵ | 4 |
| Dense continuations, five χ | 4,000 | L = 24 | 2.5×10⁴ | 16 |
| Committor relaunches, χ = 1.5 | 4,000 | L = 24 | 10³ / 3×10³ | 300 / 72 |
| Fluctuation–response twins, four χ | 4,000 | L = 24 | 10³ | 48 |
| Small-timestep re-equilibrations | 4,000 | L = 24 | 5×10² | 42 |

Continuation runs extend selected conditions to t̂ = 8×10⁵, and the replicate to 7.2×10⁵ or
1.08×10⁶. The last four rows are the kinetic, fate, response and calibration runs described below. Paper-scale runs store 200 snapshots. Analyses use the late half of each run.

Tier-membership measurements (Section 3.10) add: entropy production probes at dt = 6.25e-4 from
equilibrated configurations, with a paired same-configuration control; Onsager cross-coefficients
from 40 configurations equilibrated at the operating point (chi, chi2) = (0.5, 0.5), measured by
central differences with common random numbers at four window lengths; and effective-temperature
measurements from 200 starts per condition, with mobility from the drift under a small force on one
species and the conjugate fluctuation from the collective coordinate of that species.

The kinetic and fate measurements of Sections 3.5 and 3.9 use dense continuations (Δt = 50
between snapshots, t̂ = 2.5 × 10⁴) of the late configurations at χ ∈ {0, 0.25, 0.5, 0.75, 1, 1.5}
× 3 seeds and one monodisperse seed at N = 4000, and 300 committor launches of t̂ = 1000 plus 72
of t̂ = 3000 from the χ = 1.5 dense configurations. The single-particle fluctuation–response
measurement of Section 3.10 uses 48 twins of t̂ = 1000 at χ ∈ {0, 0.5, 1, 1.5} (random-sign
forcing of every particle at f = 10⁻⁵, four starts × 3 seeds), and the configurational
dissipation estimator uses small-timestep re-equilibrations of 18 stored runs. All of it is at
N = 4000.

The comparisons of Section 3.11 add 184 *Chlamydomonas* axonemes at 1000 frames per second across
eight ATP concentrations, 272,974 frames, from ref. [6], and a separate Vicsek flock simulation
(Section 2.6) run in four variants at six seeds each. Section 3.12 uses four published pairwise
sociomatrices and one published behavioural table, and adds no simulation.

**Convergence.** Equilibration time grows steeply with system size, and an underequilibrated
large box systematically *understates* the largest cluster — mimicking a finite characteristic
cluster size. Every continuation run raised the largest-cluster fraction: +3.4% (N=4000, χ=1.5),
+8.8% (N=9000, χ=1.5), +50.9% (N=16000, χ=1.5), +86.8% (N=9000, χ=1), and +68% on average at
N=22,000 (four seeds, range +24% to +115%). No scale-selection claim below rests on a run that has
not been continued and shown to plateau. Collapse of the within-run and between-seed scatter is a
more reliable convergence indicator here than a slope fit, because cluster counts fluctuate by up
to 30% within a run. That collapse is clean at N=4,000 and N=9,000 (between-seed spread in lcf
falling to ±0.001 at N=9,000) but not at N=16,000, where two converged seeds still differ (0.709
and 0.857); the N=16,000 point in Section 3.8 therefore carries a correspondingly wide error, which
its role in the scaling fit must be read against.

### 2.6 The comparison systems

Three systems outside the colloid model are used as controls and comparisons. Because the point of
each is that the *same* estimator is applied, we specify the shared pipeline once.

**The circulation estimator.** For any two observables (a, b) sampled at uniform intervals, the
signed area rate is A = ½⟨a ḃ − b ȧ⟩, computed from successive differences. At stationarity,
detailed balance implies A = 0 for any observable pair, so a nonzero A is sufficient evidence of
broken time-reversal symmetry; a zero A is not evidence of detailed balance, since the current may
lie outside the projection (Section 2.4). Significance is assessed against 80 phase-randomised surrogates per object, which preserve
the power spectrum of each channel while destroying cross-channel phase relations, and reported as
|z| = |A − ⟨A_surr⟩| / sd(A_surr). We take absolute values because the sign of a principal
component is arbitrary, so signed z-scores are not comparable across objects.

**Cilia.** Axoneme centrelines from ref. [6] are converted to a tangent-angle representation
ψ(s, t) along the arclength s, which removes rigid translation. Rigid *rotation* is removed by
subtracting the per-frame mean angle, ψ → ψ − ⟨ψ⟩_s; omitting this step leaves the whole-axoneme
rotation in the signal and yields beat frequencies an order of magnitude below the known 20–50 Hz,
which is how we detected the omission. The de-rotated ψ is reduced by singular value decomposition
to its two dominant shape modes, and the estimator is applied in that plane.

**Single colloidal clusters.** To match the level of description, individual clusters are tracked
over a 101-snapshot window and represented by a low-order shape descriptor: the deviatoric part of
the gyration tensor together with the normalised third and fourth mass multipoles. This is the
closest available analogue of a tangent-angle representation for an object with no centreline. The
descriptor is reduced by principal component analysis to two modes and the identical estimator and
surrogate protocol applied. Only clusters persisting for the full window are used, which limits the
comparison to three to six clusters per condition.

**Vicsek flocks.** A standard Vicsek model is run with N = 800 in a periodic box L = 20,
interaction radius r = 1, speed v = 0.3, dt = 1, angular noise amplitude η = 0.25, for 6000 steps
of which the first 2000 are discarded:

    θ_i(t+dt) = arg⟨ e^{iθ_j} ⟩_{j ∈ N_i} + ω dt + η·U(−π, π)

Four variants are used: reciprocal with the full 2π interaction neighbourhood (ω = 0); nonreciprocal
via a vision cone of half-angle 1.2 rad, so that *i* may align to *j* while *j* does not see *i*;
reciprocal with an intrinsic turning rate ω = 0.02 (chiral); and chiral with a vision cone. Six
seeds per variant. The coarse observable is the mean heading vector, whose two components form the
plane in which the estimator is applied — matched in spirit to the two shape modes used for the
axoneme — and polar order |⟨e^{iθ}⟩| is reported alongside. Turning rates
ω ∈ {0.005, 0.01, 0.02, 0.04} are additionally used to test whether the measured area rate tracks
the imposed rate.

**Dominance networks.** For a group of n individuals the net directed interaction flow
f_ij = −f_ji is an antisymmetric edge function on the complete graph. Filling all triangles leaves
no harmonic component, so the discrete Hodge decomposition is exact and two-part: f = grad(s) +
curl, with the gradient part derivable from a scalar rank potential s (obtained as s = −⟨f⟩_row)
and the curl part genuinely non-gradient. We report the cyclic fraction ‖curl‖/‖f‖. The
decomposition is applied to **log-odds**, not to raw net counts: a Bradley–Terry hierarchy has
logit(p_ij) = r_i − r_j exactly, so its log-odds flow is exactly a gradient, whereas the net-count
flow is a logistic function of the rank difference and a linear decomposition of counts therefore
retains a residual cyclic fraction that does not vanish with more data (Section 3.12).

Implementation details, since they affect the numbers. Win probabilities are formed with a
Jeffreys-type additive smoothing, p_ij = (w_ij + ½)/(w_ij + w_ji + 1), which keeps the log-odds
finite when a dyad is a clean sweep; dyads with no recorded interactions are set to zero flow and
so contribute nothing to either component. The potential is the least-squares solution
s = −⟨f⟩_row with gauge Σs = 0, and the norm is the unweighted Frobenius norm on the edge flow, so
each dyad counts equally regardless of how many interactions it carries. Matched floors are
computed from 200 regenerated Bradley–Terry replicates at the observed per-dyad counts. All four
published matrices are complete comparison graphs. Rank uncertainty is not propagated into the
floor, which is a limitation of the reported z-scores.

## 3. Results

### 3.1 Two previously proposed knobs do not vary reciprocity

The coupling-strength parameter α̂ **cancels exactly** from the antisymmetric/symmetric ratio
|a_i−a_j|/|a_i+a_j| (notation of Section 2.2): both sectors scale linearly in it. Numerically the ratio is identical to
six decimal places across α̂ ∈ {0.003, 0.005, 0.010, 0.015} at every separation tested. Sweeping
α̂ produced edge turnover *falling* from 0.601 to 0.282 — opposite to the prediction it was
meant to test — because α̂ raises overall attraction and condenses the system (n_cl 40 → 9).

The steric size-ratio sweep also fails, because the protocol of ref. [1] pins the EHD radii
while varying steric radii; it therefore varies packing, not reciprocity. Turnover falls 506×
toward the tail-large regime, but by total condensation into a single cluster (n_cl = 1.0,
lcf = 1.00), not by graded suppression of fission.

This is a caution specific to simulation: in experiment a particle's EHD and steric radii are
the same physical size, so the experimental size ratio does tune nonreciprocity.

### 3.2 Nonreciprocity produces the fragmented morphology termed arrested coarsening

Paper scale, 3 seeds per level, all other parameters fixed:

| χ | 0 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 | mono ref |
|---|---:|---:|---:|---:|---:|---:|---:|
| n_cl | 11.1 | 2.2 | 17.3 | 56.5 | 95.2 | 137.3 | 11.2 |
| lcf | 0.301 | 0.881 | 0.904 | 0.800 | 0.791 | 0.753 | 0.261 |
| Q | 0.930 | 0.911 | 0.910 | 0.909 | 0.906 | 0.901 | 0.936 |
| σ²_v | 4.19e-10 | 4.28e-10 | 3.65e-09 | 1.37e-08 | 2.76e-08 | 5.57e-08 | 3.07e-10 |

Rank correlations across χ (18 runs): n_cl ρ = +0.931; σ²_v ρ = +0.975 (the nominal p-values,
2 × 10⁻⁸ and 7 × 10⁻¹², treat the 18 runs as independent, which the blocked design makes them
not, so we do not rest on them).
χ = 0 → 1.5: n_cl 11.1 ± 0.5 → 137.3 ± 0.5 (a factor of 12.4); σ²_v × 133. We emphasise blocked
contrasts and seed-level variability because there are only three independent seed blocks and
the nominal run-level tests do not account for their dependence: with three blocked seeds per
level, and within-run fluctuations of up to 30% (Section 2.5), formal significance tests on run-level means
produce implausibly small p-values that reflect the variance-reduction of the blocked design rather
than the evidence. The effects are large enough not to need them. Both contrasts *strengthen* at paper scale relative to pilot scale, where the same ratios are
9.8× and 54×, indicating the pilot runs had not reached steady state.

A word on what "arrested coarsening" denotes here, since Sections 3.5 and 3.8 will appear to
contradict this heading. What nonreciprocity produces is the *finite-time fragmented morphology*
conventionally given that name: many clusters coexisting and continuously reorganising, on the
timescale at which the system is observed. It does so while simultaneously *unjamming* the
reciprocal gel and **accelerating** growth of the largest phase (Sections 3.5 and 4.5), and the
fragmented state is a long-lived transient rather than an asymptotic one (Section 3.8). The
reciprocal case, not the nonreciprocal one, is the genuinely arrested state. Both statements are
true and they refer to different observables — the fragment population and the majority phase.

**The antisymmetric sector therefore alters equal-time structure.** A word on terminology, which
matters because the paper also establishes kinetic trapping and long-lived transience. We use
*equal-time structure* for the observable class (n_cl, lcf — quantities computed from a single
snapshot), *late-time morphology* for what a finite-duration simulation actually realises, and
*stationary distribution* only for a genuine stationary ensemble, which at χ > 0 we do not claim to
have sampled. The result here is about the first: n_cl and lcf are single-snapshot
observables; they move by an order of magnitude under a change confined to the antisymmetric
coupling. The premise that a solenoidal coupling cannot alter the stationary density holds only
when ∇·(ρ_eq **v**) = 0, which generic nonreciprocal couplings do not satisfy. Our measurement is
of equal-time structure at late times rather than of a verified stationary density, and the two
must be kept apart: two dynamics can share a stationary distribution and approach it at very
different rates, so from a common nonequilibrium initial condition they can show very different
equal-time structure throughout any accessible window. Section 3.5 supplies exactly that
alternative — one dynamics escapes a kinetic trap that the other cannot. What the structural
contrast establishes on its own is therefore that scaling the antisymmetric coupling changes the
finite-time structure and the kinetics of condensation at fixed symmetric pair law. Whether it
also changes the stationary density is a separate, mathematical question, and it can be answered
directly for this model. The reciprocal stationary density is ρ_eq ∝ e^{−U/T}, the added drift is
**v**_a = μ**F**_a, and the drift preserves ρ_eq iff ∇·(ρ_eq **v**_a) = ρ_eq[∇·**v**_a −
**v**_a·∇U/T] vanishes identically. Writing the radial coefficients of
Section 2.4 explicitly — **F**_a,i = **F**_a,j = h(r)**r**_ij with h = −(a_i − a_j)/2 the
antisymmetric coefficient, and **F**_s,i = −**F**_s,j = h_s(r)**r**_ij with h_s = f_spring − ḡ
the symmetric one — the bracket for an isolated L–S pair is

    B(r) = (μ_L − μ_S) [ (2h + r h′) + h h_s r²/T ] ,

in units of inverse time (the same units as the per-particle dissipation table of Section 3.4).
It is not identically zero as a function of separation: B = −3.3 at r = 0.25, +0.019 at r = 0.28,
+7 × 10⁻⁴ at r = 0.5 and +3 × 10⁻⁵ at r = 0.9, changing sign where the steric spring engages at
r = s_L + s_S = 0.278. It vanishes identically only for equal mobilities, in which case the pair
propulsion is a rigid translation of the pair. Not identically zero is the relevant criterion,
and it is met. Over the
χ = 0 late configurations its per-particle value has mean +0.001 ± 0.0008 (the integration-by-parts
identity, Section 2.4) but a standard deviation of 0.006, i.e. 21 % of ã, across configurations.
The antisymmetric drift of this model does not preserve the reciprocal stationary density; it is
not a density-preserving perturbation, and the solenoidal premise fails for it as a matter of
algebra and not only of finite-time observation. The two results are kept distinct: the
density-preservation calculation establishes a property of the model's generator, while the
simulations demonstrate changes in finite-time structure and condensation kinetics and do not
establish convergence to the nonreciprocal stationary distribution.

### 3.3 Newman modularity is degenerate on these networks

Across the 14 pilot-scale conditions, cluster count spans 40× (1 → 40) while Q spans 15%, and
the two are uncorrelated: Spearman ρ = +0.31, p = 0.28. A 3-cluster state (Q=0.869) and a 33-cluster
state (Q=0.854) are indistinguishable. Louvain's resolution limit merges communities below
~√(2E) ≈ 60 nodes [4], and the degeneracy survives a 4× increase in system size. This is a known
hazard of modularity maximisation rather than a peculiarity of our graphs: the modularity landscape
is characteristically degenerate, with many structurally distinct partitions at near-identical Q
[23], which is among the reasons inferential rather than descriptive community detection is now
preferred [24].

Modularity excess over the degree-preserving null **decreases monotonically with
nonreciprocity** (0.440 at χ=0 → 0.383 at χ=1.5) and is *highest* in the monodisperse
equilibrium reference (0.462). The prediction that motivated this measurement — that constrained
organisation shows itself as a sustained modularity *excess* over a null rather than as high
absolute modularity [30] — is therefore not merely unsupported when carried over to nonreciprocal
coupling; the ordering is inverted.

### 3.4 No circulation in coarse observables, but genuine microscopic irreversibility

The signed area rate was measured in four observable planes at every χ, at both system sizes.
Every value is consistent with zero and none is distinguishable from the χ=0 equilibrium null
(all p ≥ 0.38 at paper scale, where nsnap=200 doubles the statistics).

The excess dissipation over the reciprocal control, by contrast, is clearly positive above a
threshold. The configurational estimator of Section 2.4, on re-equilibrated late configurations
at N = 4000 (three seeds per χ), gives per particle:

| χ | 0.25 | 0.5 | 0.75 | 1 | 1.5 | 2 | 3 | 5 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| excess (×10⁻³) | −0.03(2) | 0.24(3) | 1.61(17) | 2.32(22) | 7.40(32) | 12.8(3) | 23.4(16) | 64.6(34) | 178(10) |
| ÷ χ² (×10⁻³) | −0.5 | 0.95 | 2.9 | 2.3 | 3.3 | 3.2 | 2.6 | 2.6 | 2.8 |
| screening | 1.02 | 0.97 | 0.90 | 0.91 | 0.87 | 0.88 | 0.90 | 0.91 | 0.91 |

(rows: excess dissipation per particle, its ratio to χ², and the screening fraction
$-\tilde b/(\chi\tilde a)$; the χ ≥ 3 columns are production-timestep values, where the
calibration shifts the excess by less than its error). The paired trajectory probe of Section 2.4, an independent estimator, gives −0.06, −0.02,
0.80, 2.0 and 6.5 × 10⁻³ at χ = 0.25 to 1.5 at N = 1000, and 12.2 × 10⁻³ at χ = 2 and N = 4000:
consistent where it resolves, and unable to resolve the low-drive rows.

Three features of the table carry the argument of Section 3.10. The excess per χ² sits on a
plateau of 2.3–3.3 × 10⁻³ from χ = 0.75 to χ = 8; at χ = 0.5 it is 0.35 ± 0.05 of the value that
plateau would give, and at χ = 0.25 it is unresolved, −0.03 ± 0.02 × 10⁻³ against a
plateau-extrapolated 0.17 × 10⁻³. We say *unresolved below the crossover* and not zero: a
finite-noise system can carry a small activated dissipation below an apparent yielding crossover,
and the measurement does not distinguish that from an exact threshold. The unopposed term
$\tilde a$ is nearly constant in χ (0.025–0.030; it counts L–S contacts, not drive), so the
χ-dependence sits in the ratio $-\tilde b/(\chi\tilde a)$, which is 0.97–1.02 at χ ≤ 0.5 and
0.87–0.91 at every drive from 0.75 to 8: **the configurational work-rate decomposition cancels
87–91 % of the unopposed antisymmetric contribution above the crossover and, within error, all of
it below.** Which force does the cancelling changes across the crossover. Separating
$\tilde b$ into its steric, symmetric-electrohydrodynamic and divergence parts, at χ = 1.5 the
steric contact force carries 99 % of it (−0.034 of −0.034, with the symmetric EHD attraction at
−0.002 and the divergence term at +0.002), and at χ = 8 the same (107 %, −6 %, −1 %); at χ = 0.5
the symmetric EHD attraction carries 92 % (−0.013 of −0.014) and the contact force 25 %, with the
divergence term at −17 %. Above the crossover the push is held by compressed contacts; below it
the pair sits in the electrohydrodynamic well and the well holds it. At N = 1000, 9000 and 16 000
the plateau per χ² is 3.7, 2.4 and 2.4 × 10⁻³ at χ = 1.5, so the pilot-scale value is about 30 %
high and the larger boxes agree with N = 4000.

Two more cautions attach to the low-drive rows. First, what the estimator returns is the
nonconservative work rate χ⟨**F**_a∘**ẋ**⟩, and the medium heat is that plus −dU/dt; only the
total entropy production, which adds the system entropy change, is non-negative by theorem, and
away from stationarity neither the medium heat nor a difference between probes inherits that
bound automatically. We use the identity ⟨**F**_s∘**ẋ**⟩ = −d⟨U⟩/dt with −dU/dt measured at
1.5 × 10⁻⁶ per particle in units of T per unit t̂ in these windows (against table entries of order
10⁻³ in the same units), so the distinction is numerically immaterial here, but the low-drive
ensembles are aging and we do not claim a stationary identity holds for them exactly. Second, the
small residual is what has to converge with the timestep, not the large cancelling terms. On the
seed-1 configurations the χ = 0.25 excess is −0.20 (production dt), −0.03 ± 0.10 (dt = 0.005) and
+0.03 ± 0.09 × 10⁻³ (dt = 0.0025); at χ = 0.5, −0.01, +0.21 ± 0.18 and +0.40 ± 0.19; at χ = 0.75,
+0.67, +1.94 ± 0.26 and +1.33 ± 0.29, with the stationarity residual falling 4.5 → 0.34 → 0.18.
The two small-dt values agree within their errors at each χ; the statistical error of a single
seed at dt = 0.0025 (±0.1–0.3 × 10⁻³) is the larger of the two uncertainties, and the numbers
in the table above, at dt = 0.005 over three seeds, should be read with a systematic uncertainty
of that order at χ ≤ 0.75.

**The dynamics are genuinely irreversible, and that irreversibility does not appear as
circulation in any network observable measured.** Dissipation is the standard measure of how
far an active system sits from equilibrium [19, 20]; the point here is that a positive value of it
does not guarantee a signature in any given coarse observable.

### 3.5 The reciprocal reference is a kinetically arrested gel

At χ=0 the system plateaus at lcf = 0.30 from t≈1.6×10⁵ onward — it does not phase separate.
Three lines of evidence identify this as kinetic arrest rather than equilibrium:

1. Integrating the symmetric pair force gives well depths of −7.0 kT (L-L), −6.8 kT (L-S),
   −4.8 kT (S-S) against kT_eff = σ²/2. Bonds essentially never break thermally.
2. Cluster diffusion scales as D₁/n: a 1000-particle cluster moves 0.18 in a box of 24 over an
   entire run. Coalescence stalls.
3. **Decisive:** the condensed configuration evolved at χ=0 with no activity remains condensed
   (lcf 0.910 → 0.910; per-seed trajectories 0.763→0.763, 1.000→1.000, 0.966→0.967). The
   condensate is equilibrium-accessible; χ=0 simply cannot reach it. The same dispersed
   configuration is frozen at χ=0 (n_cl 10.3→10.0) but coarsens at χ=0.25 (10.3→8.7).

Weak nonreciprocity unjams this. At matched cluster size (40–120 particles) active clusters move
**10.7× faster** than thermal ones. This is an amplitude effect rather than a change of scaling,
and the evidence for that is the coherence rather than the drift. The pair propulsions add
incoherently at every size: the coherence C ≡ |Σ**p**|/Σ|**p**| falls as N^−0.53 across the full
range of cluster sizes observed, and C·√N varies only between 1.7 and 2.7 from N ≈ 5 to N ≈ 80. The
net drive therefore grows as √N and the drift speed should fall as N^−1/2 — the same exponent as
thermal diffusion, so activity cannot select a size by outrunning diffusion at one. The measured
drift exponent is consistent with that over the full size range (−0.38 per cluster, −0.47 on binned
means) but is not resolved over the small-to-medium range alone (−0.09 to −0.12), where the dynamic
range in N is too small to fit an exponent. We therefore rest the claim on the coherence
measurement, which is direct, and not on the drift fit.

Consequently cluster count is **non-monotonic**: n_cl dips to 2.2 at χ=0.25 (below the
reciprocal 11.1) as weak activity completes the arrested phase separation, then rises to 137 as
strong activity drives fission. Fluctuations peak in this crossover region (sd/mean of n_cl:
7.2% at χ=0, 29.7% at χ=0.25, 18.6% at χ=0.5, 7.1% at χ=1.5).

**The dip has a kinetic account.** Dense continuations (Δt = 50, three seeds per χ) resolve
the individual events that make and remove fragments — clusters of size ≥ 2 other than the
condensate. Tracking every cluster to its plurality successor, the condensate's *shedding rate*
(pieces of size ≥ 2 detaching per interval) and the fragments' *reabsorption probability* (per
fragment per interval, successor is the condensate) are:

| χ | 0 | 0.25 | 0.5 | 0.75 | 1 | 1.5 |
|---|---:|---:|---:|---:|---:|---:|
| fragments per snapshot, observed | 9.2 | **1.1** | 15.0 | 57 | 91 | 140 |
| shedding rate k_shed | 0.063 | 0.088 | 0.70 | 2.15 | 3.23 | 2.96 |
| reabsorption probability p_abs | 0.0073 | **0.082** | 0.047 | 0.039 | 0.038 | 0.023 |
| fragments predicted, k_shed / p_abs | 8.6 | 1.1 | 14.9 | 55 | 86 | 129 |
| fission probability per interval, S = 10–29 | — | — | 0.10 | 0.21 | 0.25 | 0.30 |

At every χ, including the reciprocal reference, the fragment count equals the shedding rate
divided by the reabsorption probability to within 8 %. We are careful about what that is. The
reabsorption probability is absorptions per fragment-interval, so k_shed/p_abs = n_f holds
exactly when the condensate sheds as many fragments as it absorbs — it is the condensate's
exchange balance, not a prediction from independent quantities — and it can only *set* n_f if the
other channels that change the fragment number net to zero. A complete budget over the same
intervals (seed 1 at each χ; `fragment_budget.py`, closure to within 0.26 events per snapshot)
shows that they do, channel by channel:

| per snapshot, χ = | 0 | 0.25 | 0.5 | 0.75 | 1 | 1.5 |
|---|---:|---:|---:|---:|---:|---:|
| shed / absorbed | 0.07 / 0.08 | 0.06 / 0.06 | 0.63 / 0.61 | 2.08 / 2.12 | 2.92 / 2.98 | 2.74 / 2.75 |
| fission / fusion | 0.17 / 0.18 | 0.03 / 0.03 | 1.27 / 1.25 | 5.45 / 5.33 | 7.58 / 7.34 | 12.1 / 12.0 |
| formation / dissolution | 0.03 / 0.02 | 0.00 / 0.00 | 1.03 / 1.03 | 6.36 / 6.43 | 12.6 / 12.9 | 19.1 / 19.3 |
| fragments | 8.3 | 2.1 | 18.5 | 61 | 98 | 138 |

Fission and fusion among fragments, and formation of dimers from monomers against complete
dissolution, are each several times larger than the condensate exchange at χ ≥ 0.75, and each
pair balances separately to within a few per cent. The source–sink accounting closes over the measured window and supports a fragment population
that is approximately stationary there, with each channel pair approximately balanced. Separate
balance of aggregate rates does not by itself establish that any one pair determines the count
— all of the rates depend on the same population and morphology — so what follows is a kinetic
interpretation rather than an independently validated predictive law. The non-monotonic fragment
count is accompanied by distinct drive dependences of condensate shedding and per-fragment
reabsorption, and the two do not switch on together:
reabsorption is already eleven times its reciprocal value at χ = 0.25, because weak activity
mobilises fragments and condensate alike, while shedding is still at its reciprocal level, so
fragments are cleared faster than they are made and the cluster count falls fivefold; shedding
then rises eightfold between χ = 0.25 and 0.5 and threefold more to 0.75, the fission probability
of an existing fragment rises from zero to 0.30 over the same interval, and the fragment count
climbs two orders of magnitude. Neither rate is monotone over the whole range — reabsorption
peaks at χ = 0.25 and declines, shedding saturates above χ = 1 — but their onsets are at different
χ, which is the kinetic account of the non-monotonic n_cl. Structure here is described by a balance
of kinetic rates, and no static functional enters the account. The budget is from one seed per χ and one window; it
should be read as demonstrating closure and the ordering of the channels, not as a precision
measurement of any of them.

We note for Section 3.10 that the shedding onset, between χ = 0.25
and 0.75, is where the excess dissipation of Section 3.4 rises to its quadratic plateau.

### 3.6 Attributing the reciprocal/nonreciprocal difference

The monodisperse-versus-bidisperse comparison confounds reciprocity with polydispersity. The
χ=0 control separates them (paper scale):

| observable | mono | χ=0 | χ=1 | share attributable to reciprocity |
|---|---:|---:|---:|---:|
| edge turnover | 0.148 | 0.260 | 0.378 | **52%** |
| σ²_v | 3.07e-10 | 4.19e-10 | 2.76e-08 | **99.6%** |

σ²_v is a clean reciprocity probe; edge turnover is about half polydispersity.

### 3.7 Observable reliability

Edge Jaccard turnover fails as a dynamical probe on four counts: about 70% of its χ=1 value
survives at χ=0 where the dynamics obey detailed balance and carry no stationary current (69% at paper scale, 71% at pilot scale); it
correlates with morphology (Spearman ρ = +0.77 at pilot scale, +0.62 at paper scale, against n_cl);
its lag-1 value is dominated by hard-cutoff flicker (0.274 at equilibrium, rising to 0.475 by
lag 32 without plateau); and it is **non-monotonic in χ at paper scale**, peaking at χ=0.75.
A monotonic trend observed at pilot scale did not survive equilibration.

σ²_v as implemented is the variance of *speed* over one snapshot interval (20,000 integration
steps, during which particles move ~0.6 radius); it is a mobility probe, not an instantaneous
velocity variance, and is not directly comparable to Fig. 2d of ref. [1].

### 3.8 No characteristic cluster size is selected; the small-cluster state is transient

The E-vs-C tradeoff account of arrested coarsening requires the system to select a finite cluster
size, which appears as an interior peak in the cluster-size distribution. Testing this requires
growing the box at fixed density: a finite characteristic size S* holds the largest cluster
constant (lcf ∝ S*/N ∝ 1/L²), whereas phase separation makes it grow ∝ N.

At χ=1.5, densities matched exactly at 6.944 particles per unit area, all points continued to
convergence:

| N | L | largest cluster | ratio | (N ratio) | lcf |
|---:|---:|---:|---:|---:|---:|
| 4,000 | 24 | 3,116 | 1.00× | 1.00× | 0.779 |
| 9,000 | 36 | 7,221 | 2.32× | 2.25× | 0.802 |
| 16,000 | 48 | 12,581 | 4.04× | 4.00× | 0.786 |

*A note on which lcf is which.* Largest-cluster fractions for χ = 1.5, N = 4000 appear at two
values in this paper and the difference is duration, not disagreement. Section 3.2's table reports
0.753, the late half of the production runs at t̂ = 4×10⁵. The value here, 0.779, is the late half
of the same runs continued to t̂ = 8×10⁵, and is the converged one; it is also the value the
rectangular-box control of Section 5 is compared against. All scale-selection claims use the
continued runs.

**Largest cluster ∝ N^1.00 ± 0.07**, where the uncertainty is the spread over all twelve
combinations of one seed per system size (range 0.93 to 1.07), against 1.00 for phase separation
and, for a finite characteristic size, the subextensive growth of the sample maximum — a
distribution with an exponential tail would give a largest cluster growing as log N, a ratio of
about 1.2 over this range of N, against the 4.04 measured. We do not treat this as a precisely resolved
exponent: it rests on three system sizes, and the two N = 16,000 seeds converge to appreciably
different largest-cluster fractions (0.709 and 0.857). The robust statement, which needs no fit, is
that the largest-cluster fraction stays approximately constant — mean 0.79, range 0.71 to 0.86
across all seeds and sizes — while N increases fourfold. That is inconsistent with an N-independent
characteristic cluster size over the range tested, and consistent with macroscopic phase
separation. The fragment
population is extensive (n_cl ∝ N^0.96) with median cluster size 3 at every box. No interior
peak appears at any size tested.

**The replicate at the source specification.** Composition follows the published
5,000:17,000 (22.7% type-I) rather than the 25% of the other series; type-I particles carry 5.1× the
EHD strength, so the higher fraction biases toward condensation:

Four seeds, each carried forward until its largest-cluster fraction stopped drifting:

| seed | t̂ = 0 – 3.6×10⁵ (**source duration**) | 3.6 – 7.2×10⁵ | 7.2×10⁵ – 1.08×10⁶ | drift in final window |
|---:|---:|---:|---:|---:|
| 1 | 0.381 | **0.811** | — | +1.1% |
| 2 | 0.479 | **0.853** | — | +0.8% |
| 3 | 0.441 | 0.536 | **0.549** | +0.0% |
| 4 | 0.554 | 0.708 | **0.870** | −0.0% |
| mean | **0.464 ± 0.037** | | **0.771 ± 0.075** | |

Seeds 3 and 4 were still drifting after the second window and were therefore extended to a third.
All four then show no detectable drift over their final observation window. We distinguish that
from ensemble or asymptotic convergence, which we cannot demonstrate: a within-run plateau in a
system whose cluster diffusion has become extremely slow is consistent with a genuine steady state
and also with arrest at a configuration-dependent macrostate. The seed heterogeneity below suggests
the latter, so we describe these as late-time plateaus rather than asymptotic values. The paired increase from the source duration to the
converged state is **+0.307 ± 0.070, a 66% rise** (paired t = +4.4, p = 0.022), and every seed
moves in the same direction, by between 24% and 113%.

One feature deserves comment: **the converged state is not unique**: three seeds settle at
0.81–0.87 while seed 3 settles at 0.549 and stays there. That heterogeneity is consistent with the
kinetic-arrest picture of Section 3.5 — a configuration that has separated into two large clusters
which then cannot find each other has no route to a single condensate on any accessible timescale.
The late-time largest-cluster fraction is therefore configuration-dependent over the windows
reached, which is itself a statement about arrest rather than about coarsening; we do not claim
the asymptotic states differ.

Two results follow. First, **composition does not explain the discrepancy**. At convergence the
22.7% replicate reaches lcf = 0.771 ± 0.075, squarely within the 0.78–0.79 seen across our 25%
runs at N = 4,000 to 16,000, rather than the 47% shortfall the underequilibrated comparison
suggested. This is a coarse comparison rather than a controlled one — the replicate differs from
those runs in system size and in χ as well as in composition — so it cannot rule composition out;
what it shows is that the corrected-composition replicate also develops a large condensate, with
a late-time fraction a few per cent from the 25 % runs and not in the direction that would
reconcile us with the source. Second, **at the
source paper's own simulation duration this reimplementation reproduces its reported state** — a
majority of particles outside the largest cluster, median cluster size 3 — and the same system
continued to twice that duration phase-separates.

The small-cluster state is therefore a long-lived transient of this model on the timescales
reached rather than its steady state — the plateaus are within-run and the seed heterogeneity
means we cannot call them configuration-dependent asymptotic states, only late-time plateaus.
This is both a validation of the reimplementation (it reaches the published state under
the published conditions) and a limitation of the model (that state does not persist).

**The transient's lifetime scales linearly with N.** Chaining each run with its continuations, the
time at which the largest-cluster fraction first reaches 0.5 is 5.5 × 10⁴, 1.44 × 10⁵ and
2.52 × 10⁵ at N = 4000, 9000 and 16 000 (χ = 1.5), and 5.7 × 10⁴ and 1.26 × 10⁵ at N = 4000 and
9000 (χ = 1): t_½ ∝ N^{1.0–1.2}. The whole lcf(t) curve collapses in t/N and not in t/√N —
at t/N = 10, 20 and 40 the three system sizes give lcf = 0.47/0.44/0.32, 0.55/0.60/0.55 and
0.77/0.77/0.75 — so the crossover time is t_x ∝ L², a diffusive crossing of the box, measured over
a fourfold range of N. The growth law behind it, from the mass-weighted mean cluster size
S_w = ΣS²/ΣS, is S_w ∝ t^{0.5–1.0} in the growth phase (t = 10⁴–2 × 10⁵) at χ = 1–1.5, saturating
once the condensate has formed; the N = 22 000 replicate shows S_w ∝ t^{0.55} over 95 snapshots
with no saturation before lcf = 0.5. At χ = 0 the same measure ages as t^{0.11} and never reaches
lcf = 0.5, and above χ = 3 the growth exponent falls again (0.43 at χ = 5, 0.20 at χ = 8, on
fewer than ten pre-saturation points, so provisional): strong drive slows coarsening in the
growth phase even though the condensate still forms.

### 3.9 The selected scale is a critical-size crossover, not a stable cluster size

Two further measurements resolve what the system selects. **Fragments** — all clusters except
the largest — have size statistics independent of system size: mean 5.02, 4.78, 5.34 particles
at N = 4,000, 9,000, 16,000, with an N-independent cutoff (p99 = 33, 35, 26) and a constant
fragment mass fraction of ~16%. A scale is therefore selected.

Its nature follows from resolving the size dynamics of individual clusters. At the production
snapshot interval this is impossible — a small cluster is typically absorbed within one interval,
so tracking reports condensate size rather than a growth increment. Dense continuations from
equilibrated configurations at Δt = 50 (three seeds, 209 748 cluster-observations of size ≥ 2
excluding the condensate) resolve it. Each cluster is followed to its plurality successor and the
interval is classified as fission (two or more pieces of size ≥ 2), evaporation (one piece plus
monomers), fusion with another cluster of size ≥ 2, absorption into the condensate, or none of
these. Evaporation must be kept apart from fission: counted together they produce a spurious crossing
of "split" and "merge" rates near S ≈ 7, whereas the fission probability itself rises
monotonically with size (0.03 at S ≈ 3, 0.24 at S ≈ 8, 0.34 at S ≈ 23, 0.54 above 100) while
fusion is flat at 0.17–0.26, so event rates alone do not select a scale. Two other measurements
do.

**The channel-resolved size drift changes sign.** Restricting to monomer and dimer exchange
— the Becker–Döring channel, with fission and fusion excluded — the mean size change per interval
is

| S | 2 | 3–4 | 5–6 | 7–9 | 10–13 | 14–19 | 20–29 | 30–49 | 50–99 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ⟨ΔS⟩ | +0.09 | −0.14 | −0.16 | −0.06 | +0.02 | +0.05 | +0.09 | +0.12 | +0.14 |
| s.e.m. | 0.003 | 0.004 | 0.009 | 0.012 | 0.017 | 0.022 | 0.032 | 0.041 | 0.055 |

negative from S = 3 to about 9 and positive above, with a zero at **S\* = 10.1** (68 % block-
bootstrap interval 9.7–10.7, defined in every resample). This is a zero of a *channel-conditioned*
drift, not of the full size dynamics: with all channels included and only absorption into the
condensate excluded, the mean drift is positive at every size (+1.6 at S = 2 to +34 at S ≈ 40),
because rare fusions with larger fragments dominate the mean, and with absorption included it
returns condensate size at every S, which is the artefact the coarse-interval attempt had
produced. Absorption itself is a small flat hazard of 1.5–3.5 % per interval at every size.
Dimers are the one exception to the sign pattern and gain on net.

**The committor of the full size dynamics is one half at $S_0$ ≈ 8, and the value is converged in
the horizon.** From the same configurations, relaunches with fresh noise followed every fragment
of size ≥ 3 to the first of: reaching S ≥ 20, falling to S ≤ 2, or being absorbed by the
condensate, with every size-change channel operating. Two sets were run: 300 relaunches of
t̂ = 1000 (25 noise realisations × 4 starts × 3 seeds, 23 800 fragments) and 72 of t̂ = 3000
(12 × 2 × 3, 5 664 fragments). The probability of the condensed outcome, q($S_0$) = P(grow $\cup$
absorbed) among decided fates, from the longer set:

| $S_0$ | 3 | 4 | 5 | 6 | 7 | 8–9 | 10–11 | 12–14 | 15–19 | 20–29 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| q | 0.14 | 0.25 | 0.26 | 0.21 | 0.36 | 0.58 | 0.53 | 0.63 | 0.87 | 0.97 |
| decided fates | 1899 | 939 | 594 | 360 | 309 | 332 | 269 | 252 | 169 | 180 |

with binomial errors of 0.01–0.03 per entry. The fate estimate is conditioned on a decision within
the horizon, and that matters near the crossover: shrink fates take longer than growth (median
first-passage 950 against 700 at $S_0$ = 8–11), so short horizons bias q upward and the crossing
downward. The crossing is therefore reported as a function of the horizon τ:

| τ | 250 | 500 | 1000 | 1500 | 2000 | 2500 | 3000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| $S_0$ at q = ½ | 7.3 | 7.5 | 7.8 | 7.9 | 7.9 | 7.9 | 7.9 |
| $S_0$ at q = ½, absorption excluded | 7.5 | 7.7 | 8.2 | 8.3 | 8.3 | 8.4 | 8.4 |
| fraction of $S_0$ = 8–11 fates decided | 0.20 | 0.32 | 0.56 | 0.73 | 0.83 | 0.89 | 0.93 |

It moves by 0.7 between τ = 250 and 1500 and by less than 0.05 thereafter, with 93 % of the fates
near the crossover decided at τ = 3000. The converged value is **$S_0$ = 7.9** (68 % bootstrap over
launches 7.8–8.1; over the six parent configurations 7.6–8.3), or 8.4 (8.2–9.3) when absorption
into the condensate is not counted as the condensed outcome; the 300-launch set at τ = 1000 gave
8.0 and 8.8, within 0.1 and 0.4 of these. Errors that treat fragments as independent understate
the true uncertainty, since fragments in one launch share a parent configuration and a noise
stream; the launch- and parent-level bootstraps above are the ones to quote.

The two criteria therefore bracket the scale rather than coincide: the committor of the full
dynamics is one half at $S_0$ ≈ 8 (8.4 without absorption), and the drift of the monomer-exchange
channel changes sign at S ≈ 10. We call S\* ≈ 8–10 a **critical-size crossover**. The committor
criterion, which is the one that involves every channel, is met at converged horizon, and on that
evidence the crossover is a critical nucleus of the fragment-size dynamics; the mean-drift
criterion of classical nucleation theory is met only within the monomer-exchange channel, and the
mean drift of the full process has no zero because rare fusions dominate it. Absorption into an
existing condensate is also a different event from a fragment growing into a supercritical
nucleus of its own, which is why both definitions of q are given. The crossover is a property of
the contact energetics rather than of the drive: the channel-resolved drift zero is 10.7, 9.3 and
10.1 at χ = 0.75, 1 and 1.5 (7.3 at χ = 0.5, with a tenth of the events), and Section 3.8 has
already shown the fragment scale independent of N. The fragment population (mean ≈ 5) sits below
it and, on the committor evidence, is subcritical.

**The fragment scale is independent of density as well as of system size.** Holding N = 4000 and
χ = 1.5 fixed while expanding the box dilutes the suspension from the baseline 6.944 particles per
unit area down to 1.736:

| density | 6.944 (L=24) | 4.444 (L=30) | 2.770 (L=38) | 1.736 (L=48) |
|---|---:|---:|---:|---:|
| largest-cluster fraction | 0.779 | 0.534 | 0.130 | 0.038 |
| number of clusters | 141 | 277 | 429 | 498 |
| **mean fragment size** | **5.0** | **5.6** | **6.9** | **6.3** |

The condensate is strongly density-dependent and essentially gone by the lowest density: there is a
threshold below which the majority phase does not form. We restrict the fragment claim to the three
densities at which a condensate exists, since at ρ = 1.736 there is no majority phase and
"fragments" and "all clusters" are the same population — not the regime the comparison is meant to
test. Over those three the mean fragment size moves from 5.0 to 6.9, a 38% swing across a 2.5-fold
change in density. That is a good deal weaker than "independent", and we claim only that the
fragment scale stays within a narrow band — the same band it occupies over a fourfold change in
system size (Section 3.8) — while the condensate fraction changes by a factor of six. The lowest-
density point was not continued and is the least converged of the four.

The state is thus a condensate coexisting with a population of small clusters, subcritical on the
committor evidence, that form by shedding from the condensate and by nucleation from monomers and
are removed by reabsorption and dissolution, with the crossover at S\* ≈ 8–10 separating the two
populations. No mechanism caps cluster size, consistent with §3.8: above S\* the committor rises
to one, the largest cluster grows as N^1.00 ± 0.07, and the system phase-separates.

### 3.10 The limits of three commonly used nonequilibrium diagnostics

The preceding sections establish what nonreciprocity does to this system. We now ask what kind of
variational description it admits. The question is easy to render vacuous — a framework asserting
that "some variational principle applies" forbids nothing — so we treat tier membership as an
experimental matter. We stress that the following is a working taxonomy of *commonly invoked
diagnostics*, not an established classification: the existence of a single effective temperature in
particular is a nonequilibrium construct in the sense of Cugliandolo and Kurchan [11], and is
neither necessary nor sufficient for near-equilibrium response. With that caveat:

| tier | regime | variational object | signature |
|---|---|---|---|
| I | equilibrium | free energy; MaxEnt | EPR = 0; detailed balance; Boltzmann |
| II | linear response | Rayleighian R = Φ̇ + Ψ, Ψ quadratic in rates [21] | fluxes linear in the forces; Onsager reciprocity among conjugate pairs |
| III | far from equilibrium | no general principle | nonlinear flux–force relations; reciprocity fails |

It is worth noting at the outset that least action does not fail for nonreciprocal systems in the
way sometimes supposed. The Onsager–Machlup action S[x] = ∫ (ẋ − μF)²/4D dt is well defined for
*any* drift field, gradient or not. What fails is the equilibrium corollary — detailed balance, and with it structure selected by
minimising an energy; a Boltzmann-shaped stationary density is not excluded in general, though
Section 3.2 shows this model's antisymmetric drift does not preserve the reciprocal one. That
corollary, and not least action itself, is what the solenoidal proposal of Section 1 relied upon.

We tested three tier-II signatures.

**Quadratic dissipation.** Section 3.4 gives the excess dissipation per particle at N = 4000
from χ = 0.25 to 8, on a single system size and with a single estimator. Above χ = 0.75 it is
quadratic: per χ² it is 2.3–3.3 × 10⁻³ across an elevenfold range of drive and a hundredfold
range of signal, and at χ = 2–8 alone it is flat to 4.8 %. Below χ = 0.75 it is not a quadratic
with a correction; it falls below the quadratic plateau, to a third of it at χ = 0.5, and is
unresolved at χ = 0.25 (−0.03 ± 0.02 × 10⁻³ on three seeds against a plateau-extrapolated 0.17).

We are careful about what this does and does not test, and two things have to be separated.

*This row is not a fixed-configuration measurement.* Each χ is probed from a configuration
equilibrated at that same χ, so structure co-varies with drive along the row. That makes the clean
quadratic above the threshold less trivial than a frozen-structure result would be, but also
harder to attribute, since two things change together.

*At genuinely frozen configuration the quadratic form is largely fixed by construction.* The
antisymmetric force is exactly linear in χ; the dynamics are overdamped with
configuration-independent mobility; and the symmetric sector is conservative, so its own
contribution to the Stratonovich heat is an exact differential which averages to zero in a
stationary state. What remains is $\chi(\tilde b+\chi\tilde a)$ in the notation of Section 2.4, and the quadratic
term must dominate at large χ whatever the system does. The only informative quantity is therefore
$\tilde b$, the part of the drive cancelled by the contact response.

*The construction's cross-term is contact screening, and it is quadratic, not linear.* Fitted
as aχ² + bχ over the low-drive rows, the data return a negative b; but no such form can hold at
the origin, since the excess is non-negative for every χ and zero at χ = 0, so its slope there
cannot be negative, and the fitted form would go negative below χ ≈ 0.45. The configurational
decomposition shows what such a fit describes. The unopposed
term $\tilde a$ is nearly constant in χ; $\tilde b$ is proportional to $-\chi\tilde a$ with coefficient 0.87–0.91 at every drive
from χ = 0.75 to 8, and with coefficient 0.97–1.02 at χ ≤ 0.5. An apparent linear term is a quadratic
structural screening: $\tilde b$ vanishes identically in a Boltzmann ensemble, and its measured value at
χ = 0 is indeed zero within noise (+0.0003 on $\tilde a$ = 0.030), so it is a response of the structure to
the drive and not a property of the reference. Above the crossover the configurational decomposition cancels 87–91 % of the unopposed
antisymmetric contribution, almost entirely through the steric contact force; below it the
cancellation is complete within error and is carried by the symmetric electrohydrodynamic
attraction instead (Section 3.4).

*The crossover is structural.* The cancellation steps from ≈ 1.0 to 0.90 between χ = 0.25
and 0.75, which is where Section 3.5 puts the onset of condensate shedding and fragment fission,
and where the cancelling force changes from the electrohydrodynamic well to the steric contact.
Below it the push is held by the pair well and the excess is unresolved; above it particles are
pushed onto their contacts, the condensate sheds, and the residue the contacts do not hold is
what dissipates. The
measurement is therefore best read as establishing that dissipation and fragmentation share one
threshold, and as a validation of the configurational estimator against the trajectory probe
across a hundredfold range of signal, rather than as an independent test of linear response:
above the crossover the quadratic form is built into the construction, and below it the excess
is too small to fit a form to.

**Onsager reciprocity.** This requires a second, independent drive. A uniform field applied to one
species is unsuitable: the system is isotropic, so the species drift vanishes at zero field for
every χ and the cross-coefficient is zero by symmetry rather than by physics, rendering the test
vacuous. We instead opened a second nonreciprocal channel in the collision force, with its own
parameter χ₂ and the same pair-co-propulsion structure as the electrohydrodynamic one. Entropy
production is then bilinear in two internal thermodynamic forces,

    T · EPR = χ J₁ + χ₂ J₂,    J₁ = Σ f₁ · v,   J₂ = Σ f₂ · v,

identifying (χ, χ₂) as forces and (J₁, J₂) as their conjugate fluxes. Measured at the operating
point (χ, χ₂) = (0.5, 0.5) from configurations equilibrated at that point:

| window T | symmetric part (L₁₂+L₂₁)/2 | t | antisymmetric part (L₁₂−L₂₁)/2 | t |
|---:|---:|---:|---:|---:|
| 200 | −1.08e−8 ± 6.5e−7 | −0.02 | +9.98e−6 ± 6.3e−7 | +15.8 |
| 400 | −6.81e−8 ± 8.6e−7 | −0.08 | +9.40e−6 ± 8.2e−7 | +11.4 |
| 800 | −3.96e−7 ± 9.9e−7 | −0.40 | +9.48e−6 ± 9.6e−7 | +9.9 |

The symmetric part of the cross-coupling vanishes at every window length, while the antisymmetric
part is large and stable: **the measured differential cross-response is predominantly
antisymmetric**, L₁₂ ≈ −L₂₁.

**A structurally dissimilar second channel does not reproduce this.** Because both channels above
were built with the same pair-co-propulsion structure, the antisymmetry may follow from that
construction symmetry rather than from any physical relation. We therefore repeated the
measurement with a second channel made structurally unlike the first — cubic rather than linear
size contrast, weighted toward the outer part of the overlap rather than uniform across it — from
configurations equilibrated at the operating point with that channel:

| window T | symmetric part | t | antisymmetric part | t |
|---:|---:|---:|---:|---:|
| 200 | +8.67e−7 ± 4.3e−7 | +2.04 | +2.16e−6 ± 4.2e−7 | +5.2 |
| 400 | +1.62e−6 ± 5.2e−7 | +3.12 | +1.17e−6 ± 5.1e−7 | +2.3 |
| 800 | +2.70e−6 ± 5.9e−7 | +4.57 | **+1.18e−8 ± 5.9e−7** | **+0.0** |

The pattern inverts. With dissimilar channels the antisymmetric part *decays to zero* while a
symmetric part grows and becomes significant — the opposite of the identical-channel case, where
the antisymmetric part plateaued and the symmetric part remained at zero throughout. The result
is reproduced at two equilibration lengths (T = 2×10⁴ and 1×10⁵), and neither component has
plateaued at the longest window, so the dissimilar case is not itself converged; but its direction
is consistent and unambiguous.

**We therefore do not claim antisymmetry as a general property of the cross-response between
nonreciprocal channels.** It is observed robustly when the two channels share a structure and is
absent when they do not, which is what the construction-symmetry explanation predicts. The
identical-channel measurement stands as a measurement; its interpretation as a physical reciprocity
relation does not survive this control.

We deliberately stop short of calling either measurement Onsager–Casimir reciprocity [7, 8]. That
identification would require a defined microscopic time-reversal operation with the parities of all
time-odd control parameters specified, and it cannot be inferred retrospectively from an observed
sign. Three further caveats apply independently of the control just described. The coefficients are
measured about a nonequilibrium operating point rather than about equilibrium; χ and χ₂ are
internal coupling constants in the equations of motion, not thermodynamic affinities, so the
bilinear identity above is a definition rather than a derivation; and the two-channel construction
is ours, not a feature of the source model.

The identical-channel antisymmetry admits a physical reading, and the dissimilar-channel
control was built to test it. Entropy production is the quadratic form EPR = Σ L_ij X_i X_j, to which only the
symmetric part of L contributes, so a purely antisymmetric off-diagonal block dissipates nothing;
an antisymmetric cross-coupling would therefore be *reactive* rather than dissipative, and that
would explain why the entropy production reduces to a single-coefficient quadratic. A naive parity
argument predicts the opposite — both fluxes have the form Σ f · v with f even and v odd under time
reversal, making both odd and their signatures equal — so an antisymmetry, if physical, would
demand a mechanism, for which the most economical candidate is a *gyroscopic* coupling in the
space of collective coordinates, Magnus- or Coriolis-like, entering the response matrix
antisymmetrically without contributing to dissipation. The dissimilar-channel control does not
support that reading, and we report the gyroscopic coupling as the hypothesis the control tested
and did not confirm. **A defensible Onsager test in this model would require an
equilibrium Green–Kubo derivation or an explicit Fokker–Planck computation of L_ij, and we do not
have one.**

**Effective temperature.** At equilibrium the fluctuation–dissipation theorem [25] fixes
D/μ = kT for every degree of freedom; here mobility is 1/s_i and the noise variance σ²/s_i, so
D/μ = σ²/2 exactly for both species. Out of equilibrium the ratio defines an effective temperature,
and a single T_eff shared across degrees of freedom is a tier-II signature. The collective species coordinates X = Σ_{i∈species} x_i, with the mobility from the drift under
a small force on one species, are the natural first choice, and we report that measurement below
for the centre-of-mass drift it exposes; but they cannot answer the single-temperature question,
because once the whole-system drift is removed the two species coordinates are one degree of
freedom up to sign
(X_rel(L) + X_rel(S) = 0), their apparent temperatures must differ by exactly
(n_S/n_L)(μ_S/μ_L) = 6.85 (measured 6.7), and a force applied to a whole species is not conjugate
to any single particle's fluctuation.

The measurement that does answer it is on single-particle coordinates. Every particle receives a
force f ε_i x̂ with ε_i = ±1 independent, and the run shares its noise stream with an unperturbed
twin. The paired response ⟨ε_i Δx_i⟩/f is the mean *diagonal* displacement response, cross terms
averaging out over ε, and it is conjugate to the mean single-particle mean-squared displacement;
T_eff(t) = MSD(t)/2χ(t) then equals T at every t at equilibrium. Linear response holds to t ≥ 200
at f = 10⁻⁵ against f = 10⁻⁶. Twelve twins per χ (three seeds × four late starts, N = 4000,
t̂ = 1000):

| χ | T_eff/T, t = 1–10 | t where T_eff/T > 1.1 / 1.2 | T_eff/T, t = 100 / 400 / 800 | L ÷ S, t = 400 / 800 |
|---|---:|---:|---:|---:|
| 0 | 1.00–1.04 | never / never | 1.00 / 0.98 / 0.99 | 1.00 / 1.02 |
| 0.5 | 1.02–1.04 | 216 / 340 | 1.06 / 1.28 / 1.74 | 1.01 / 0.98 |
| 1 | 1.01–1.05 | 37 / 74 | 1.33 / 2.60 / 3.94 | 1.03 / 0.93 |
| 1.5 | 1.05–1.08 | 13 / 30 | 1.81 / 4.38 / 7.17 | 1.02 / 0.89 |

The estimator calibrates exactly: at χ = 0 the ratio is unity within 3 % at every lag from 1 to
800 on both species. Out of equilibrium three things are established. First, a single effective
temperature exists and equals T over a window that shrinks with drive: FDT holds to within 10 % up
to t ≈ 220 at χ = 0.5, 37 at χ = 1 and 13 at χ = 1.5, with t_FDT falling roughly as χ⁻². This is an
operational short-lag agreement window in the time domain — a displacement-based ratio is an
integrated quantity, and agreement over short lags does not establish FDT at every frequency
above 1/t_FDT; a spectral statement would need the correlation and dissipative response spectra,
which we did not compute. Second, at longer lags the ratio grows as
(T_eff/T − 1) ∝ t^{1.0} at χ = 1 and 1.5 without a plateau to t = 800; the t^0.89 growth of the
species-coordinate measurement is this regime. Third, **the two species share the same T_eff(t)
at every lag**, to 2 % at t = 400 and 10 % at t = 800, although their bare mobilities differ by a
factor of 1.5. The single-temperature question therefore has an answer on independent
coordinates: one temperature, shared, equal to the bath temperature over a short-lag window, and
a shared function of the observation time rather than a number beyond it. This is agreement of
the one fluctuation–response relation we tested, on single-particle coordinates; it is not a
demonstration of near-equilibrium response for every observable. The violation appears
only at lags longer than t_FDT(χ), which places the dissipation of Section 3.4 in cluster-scale
motion rather than at contact scale; we did not evaluate the Harada–Sasa sum rule that would make
that quantitative.

*The species-coordinate measurement and the centre-of-mass drift.* Mobility from the drift under a small force
on one species, fluctuation from the collective coordinate, perturbed and unperturbed runs sharing
a noise stream, verified linear down to f = 3 × 10⁻⁶; at equilibrium T_eff/kT = 0.94, 1.04 and
1.09 at windows T = 200, 400 and 800.

| | χ = 0 | χ = 1.5 |
|---|---:|---:|
| system centre of mass, ⟨ΔR²⟩ | t^1.01 | **t^1.96** |
| species coordinate, ⟨ΔX²⟩ | t^1.08 | t^1.86 |
| species coordinate, drift removed | t^0.59 | t^1.47 |
| T_eff/kT at T = 200 / 400 / 800 | 0.94 / 1.04 / 1.09 | **5.18 / 9.54 / 17.66** |

About half of the growth on that coordinate is the aggregate translating under a net internal
force — a finite-size, window-limited effect rather than a new phenomenon (its amplitude falls
roughly as N^−1/2 and it decorrelates over long runs) — and the residue is internal. The
single-particle measurement above is free of that confound quantitatively rather than by
construction: a random-sign pattern still carries a net force of order f√N (|Σε_i|/√N = 0.98 in
our draws), but the centre-of-mass mean-squared displacement is below 1 % of the single-particle
value at every lag and every χ in the same runs (0.2–1.0 %), so neither the spontaneous
centre-of-mass motion nor the residual net force contributes at the level of the numbers above.
The large–small difference on species coordinates is fixed by the identity above and carries no
information about a shared temperature.

**Interpretation: tier membership is timescale-dependent, and the boundary is measured.** We set
out expecting these three diagnostics to place the system in a tier. Two of them cannot, for
reasons of construction, and the third places it in two regimes at once, separated by a lag.

The quadratic dissipation is real and cleanly measured, but weakly diagnostic on its own: above
the crossover the quadratic form is built into the construction, and below it the excess is
unresolved. What it establishes is a crossover in χ shared with the structural transition, and a
validation of the configurational estimator. The Onsager measurement returned a stable,
plateau-tested antisymmetry that the dissimilar-channel control shows to be most parsimoniously a
property of how we built the second channel; it does not survive as a physical reciprocity
relation. The effective temperature, measured on independent coordinates, is the one that returns
a definite answer: a single effective temperature equal to T exists, shared by both species, over
lags shorter than t_FDT(χ), and a shared, growing violation beyond.

The honest summary is that **the near-equilibrium description applies over a lag window that
shrinks as the drive grows, and fails beyond it**, with the failure carried by cluster-scale
motion; the agreement is established in the time domain, over lags shorter than t_FDT(χ), and we
do not claim it frequency by frequency. The position the measurements support is this: FDT with a shared temperature at contact and intra-cluster
timescales, and a growing, shared violation at cluster-motion timescales, with the boundary
between them a measured function of χ.

We none the less think the exercise was worth conducting, and record its outcome in this form
deliberately. A framework asserting that "some variational principle applies" forbids nothing; what
makes tier claims falsifiable is that each tier carries experimental signatures and that membership
be *measured*. What this section demonstrates is that measuring it honestly is harder than the
taxonomy suggests. Two of the three signatures proved near-tautological or construction-dependent
once examined closely, and the third gave opposite answers on two coordinates until the conjugate
pair was chosen correctly. That is a result about the diagnostics as much as about the system, and
it is the reason we present the taxonomy above as a working list of commonly invoked tests rather
than as a classification.

One boundary should be stated explicitly in any case, because our own data mark it. Entropy
production is quadratic *at fixed structure*, but the structural response to χ is strongly
non-linear: cluster count is non-monotonic, dipping at χ = 0.25 before rising by an order of
magnitude (Section 3.5), and Section 3.5 gives the mechanism — a balance of two kinetic rates
with different thresholds. Whatever variational description ultimately applies to this system can
at best govern its dynamics at fixed structure; *structure selection* here is a rate-balance fixed
point, and we know of no variational principle whose extremum reproduces it. Identifying that gap
seems to us more useful than obscuring it.

### 3.11 A biological comparison: does irreversibility survive coarse-graining?

The colloid suspension is useful here for a reason independent of Section 3.10's inconclusive tier
question: it is nonreciprocal, dissipative and strongly self-organising, and it is uncontroversially
not alive. Any proposed signature of biological organisation must therefore distinguish itself from
this system, not merely from equilibrium. What the comparison requires is only that the colloid be
driven, far from equilibrium and demonstrably irreversible at the microscopic level, all of which
Section 3.4 establishes independently.

Section 3.4 supplies the natural axis of comparison. The colloid suspension has provably positive
entropy production, yet the signed area rate in every coarse observable plane is indistinguishable
from the equilibrium null at both system sizes. Its irreversibility is real but microscopic: it
does not survive coarse-graining.

That such circulation exists in beating axonemes is not itself new: Battle et al. [12] demonstrated
probability-current loops in the phase space of beating *Chlamydomonas* flagella and isolated
mammalian axonemes a decade ago, establishing broken detailed balance at mesoscopic scales in
active biological systems, and the subsequent literature has developed both the measurement and its
coarse-graining caveats extensively [13, 14]. Our contribution here is not the observation of
circulation in cilia but the *controlled comparison*: the same estimator, the same surrogate
protocol, and — following the control below — the same level of description, applied to a synthetic
system whose entropy production is independently known to be positive and whose reciprocal limit is
a reciprocal reference whose stationary current vanishes.

We asked whether a biological system differs on that specific axis, using published high-speed
recordings of reactivated *Chlamydomonas* axonemes [6, 22]: 184 usable axonemes imaged at 1000 frames
per second across eight ATP concentrations from 50 to 1000 µM, 272,974 frames in total. Axoneme
shapes were converted to tangent-angle representations, removing rigid translation and rotation,
and reduced by principal component analysis to two dominant shape modes; the signed area rate was
then computed in that plane against phase-randomised surrogates, using the identical estimator
applied to the colloid observables.

We do not test whether entropy production is quadratic in the ATP drive, because that question
is ill-posed for this dataset. The chemical potential of ATP hydrolysis depends on the
concentrations of ADP and inorganic phosphate, on temperature and ionic conditions, and on
standard-state corrections, none of which are reported for these reactivation buffers; it is
therefore not fixed by the ATP concentration alone. What can be said without those quantities is
that reactivated axonemes with an ATP-regeneration system are chemically driven far from
equilibrium under all conditions in the series, so a linear-response question posed in terms of the
ATP drive answers itself rather than being decided by measurement.

| [ATP] µM | n | beat frequency (Hz) | median \|z\| | fraction \|z\|>2 | \|area per cycle\| |
|---:|---:|---:|---:|---:|---:|
| 50 | 10 | 14.5 | 3.3 | 1.00 | 6.07 |
| 66 | 15 | 23.7 | 3.3 | 1.00 | 6.24 |
| 100 | 11 | 27.7 | 3.0 | 0.91 | 6.23 |
| 240 | 19 | 46.3 | 2.7 | 0.89 | 6.15 |
| 370 | 29 | 58.4 | 3.0 | 0.97 | 5.90 |
| 500 | 13 | 91.7 | 4.5 | 0.92 | 5.95 |
| 750 | 44 | 94.8 | 3.2 | 0.89 | 5.71 |
| 1000 | 43 | 65.0 | 3.3 | 0.91 | 6.09 |
| **total** | **184** | | | | |

Pooled over all 184 axonemes, the median |z| is 3.2, with 92% exceeding |z| = 2 and 55% exceeding
|z| = 3.

**A matched-level control.** The colloid measurement of Section 3.4 uses system-averaged network
observables over ~10⁴ particles, whereas the cilia measurement uses the two-mode shape space of a
single axoneme. These are not commensurate levels of description, and the colloid null could in
principle arise from ensemble averaging alone: in an isotropic, statistically homogeneous
suspension most global scalar pairs have vanishing signed area by symmetry, and phase-incoherent
circulation across many clusters averages to zero — precisely the incoherence documented in
Section 3.5, where C√N is constant.

We therefore repeated the measurement at a matched level of description: a single tracked cluster,
represented by a low-order shape descriptor (deviatoric gyration tensor plus normalised third and
fourth mass multipoles, the closest available analogue of the axoneme tangent-angle
representation), reduced by PCA to two modes, with the identical estimator and surrogate protocol.

| system | level of description | median \|z\| | fraction \|z\|>2 |
|---|---|---:|---:|
| colloid, χ = 0 | single tracked cluster | 0.81 | 0.00 |
| colloid, χ = 1 | single tracked cluster | 0.85 | 0.00 |
| colloid, χ = 1.5 | single tracked cluster | 0.94 | 0.00 |
| *Chlamydomonas* axoneme | single organelle | **3.2** | **0.92** |

No circulation is detectable in the tracked clusters, so the contrast does not stem from the choice
of coarse-graining level. Two limitations bound this. Only three to six clusters per condition
persist long enough to be tracked over the 101-frame window, so the colloid side rests on few
objects — uniformly null, but few. And the persistence criterion may itself select against the
signal: a cluster undergoing strong shape dynamics is less likely to survive intact for the full
window, so the tracked sample may be biased toward unusually stable, non-cycling objects. The
defensible claim is that circulation is absent from the clusters we could track, not that single
active colloidal clusters never circulate.

**A synthetic positive control, and it changes the interpretation.** The comparison so far
confounds two contrasts: biological versus synthetic, and *collectively ordered* versus spatially
distributed and uncoordinated. To separate them we applied the identical estimator to a Vicsek
flock, which is synthetic and collectively ordered, in four variants:

| variant | polar order | median \|z\| | fraction \|z\|>2 |
|---|---:|---:|---:|
| standard (reciprocal, no chirality) | 0.709 | 0.55 | 0.00 |
| vision cone (**nonreciprocal**, no chirality) | 0.875 | 0.32 | 0.00 |
| chiral (reciprocal, limit cycle) | 0.015 | 1.80 | 0.33 |
| chiral + vision cone | 0.069 | **5.58** | **1.00** |

**A purely synthetic system produces circulation stronger than the cilia signal.** Neither
collective order alone nor nonreciprocity alone suffices — both give a null — but the chiral
variants circulate strongly. Note what the polar-order column says about *why*: the two chiral
variants have polar order 0.015 and 0.069, an order of magnitude *below* the two non-chiral ones.
The strongest circulation therefore occurs where the mean heading vector is smallest. These are not
strongly polar-ordered flocks, and we should not describe them as such: they are Vicsek-type
ensembles with an imposed turning rate, in which the residual mean heading rotates at that rate.
The control accordingly shows that circulation does not require biology; it does *not* show that
collective order plus chirality generates the effect.

Varying the turning rate identifies the origin: the measured area rate is 0.0056, 0.0089, 0.0141
and 0.0352 for ω = 0.005, 0.01, 0.02 and 0.04, a ratio to ω of ≈0.9 throughout. The signal is
therefore largely the *rotation of the mean heading vector at the imposed rate*, not an emergent
property.

The chiral control therefore rules out the reading that coarse-grained circulation distinguishes
biological from non-biological organisation. What the measurement detects is whether the system possesses a
**cyclic collective mode in the chosen projection** — a structural fact about the observable, which
may be emergent (the ciliary beat, arising from motor coordination) or imposed (a single-particle
turning rate), and which the estimator cannot distinguish. The colloid suspension lacks such a mode
at every level of description we examined; cilia possess one; and a synthetic flock can be given
one by construction.

What survives is the narrower and better-supported statement: the colloid suspension has provably
positive entropy production and yet no circulation in any coarse observable, at either the
system-averaged or single-cluster level. Microscopic irreversibility need not project onto
macroscopic observables. That is a statement about coarse-graining, not about life.

A second and unanticipated result emerges from the same analysis. The area enclosed per beat cycle
is 5.71 to 6.24 in standardised shape coordinates at every ATP concentration — flat across a
twentyfold range of concentration and a sixfold range of beat frequency — and that value is
approximately 2π. We are careful about how surprising this is: if two PCA modes of a periodic
waveform are standardised to unit variance and are close to phase quadrature with near-sinusoidal
profiles, an enclosed area of 2π follows almost by construction, so the value itself largely
restates a property of principal component analysis applied to periodic data. The non-trivial
content is the *invariance* — the tightness of the range, 5.71 to 6.24, across a twentyfold
concentration range — which says that the beat remains equally close to a clean quadrature limit
cycle at every drive. On that reading the beat traces a limit cycle whose *shape* is saturated and
ATP-independent, while
ATP sets only the *rate* at which the fixed cycle is traversed, with beat frequency rising from
14.5 to approximately 95 Hz in the Michaelis–Menten manner known for axonemal dynein. Geometry and
kinetics separate cleanly.

This comparison also reframes the circulation measure itself. As a probe of nonreciprocity it is
poor: our strongly nonreciprocal colloid registers none. As a probe of *biological* organisation it
is equally poor, since an imposed turning rate reproduces and exceeds the biological signal. What
it does report is whether a cyclic mode is present in the chosen projection — a structural fact
about the observable, which may be emergent or imposed, and which the estimator cannot tell
apart.

We are careful about the scope of this claim. Three systems establish that the axis exists and that
systems differ along it; they do not establish what governs position on it beyond the presence of a
cyclic collective mode, and one synthetic system with such a mode is not a survey. Circulation
measured in a two-mode projection is moreover a lower bound on irreversibility, and no system's
full phase space was searched.

### 3.12 A third system where the decomposition is exact: social dominance

Every result above probes the gradient/circulating question indirectly, through the consequences of
a decomposition rather than the decomposition itself. In the colloid the antisymmetric sector is
something we impose; whether the resulting *probability current* circulates has to be inferred from
observables. There is a class of system in which the analogous decomposition can be computed
exactly and in finite dimension, and it is worth doing because it tests the same proposition on
data we did not generate.

For a group of n animals, the net directed interaction flow f_ij = −f_ji is an antisymmetric edge
function on a complete graph. Filling the triangles leaves no harmonic component, so the discrete
Hodge decomposition [28] is exact and strictly two-part: a gradient component that *is* derivable from a
scalar rank potential — the social analogue of a conservative force — and a curl component that is
genuinely non-gradient, the rock–paper–scissors structure. The cyclic fraction ‖curl‖/‖f‖ is
therefore a direct measurement of the question the colloid work can only reach through consequences:
**is this antisymmetric coupling a gradient, or does it circulate?**

**A methodological result comes first, because it changes what the answer is.** The decomposition
must be applied to log-odds, not to raw net counts. A Bradley–Terry hierarchy has
logit(p_ij) = r_i − r_j exactly, so its log-odds flow is exactly a gradient, whereas the net-count
flow is a *logistic* function of the rank difference; a linear decomposition of counts therefore
retains a cyclic component that does not vanish with more data. That component is not
mathematically spurious — it is genuinely present in the observed count-flow field — but it is
spurious *as evidence of intransitivity*, which is the only thing anyone wants to read it as. On
synthetic perfectly transitive groups:

| interactions per pair | cyclic fraction, raw counts | cyclic fraction, log-odds |
|---:|---:|---:|
| 5 | 0.262 ± 0.093 | 0.288 ± 0.114 |
| 20 | 0.198 ± 0.054 | 0.197 ± 0.079 |
| 100 | **0.157 ± 0.030** | 0.094 ± 0.042 |
| 500 | **0.158 ± 0.013** | **0.044 ± 0.020** |

A perfectly transitive group reads as 16% intransitive on raw counts however much data is
collected. Published claims of intransitive dominance resting on a linear decomposition of counts
should be checked against this floor.

**Applied to four published sociomatrices** (bundled with the `compete` package [29]), each against a
null obtained by fitting that matrix's own rank potential, regenerating a perfectly transitive
Bradley–Terry truth at its own per-pair counts, and re-decomposing:

| dataset | n | interactions per ordered pair | cyclic fraction | matched transitive floor | z |
|---|---:|---:|---:|---:|---:|
| mouse | 12 | 3.5 | 0.495 | 0.580 ± 0.054 | −1.6 |
| caribou | 20 | 4.3 | 0.607 | 0.579 ± 0.028 | +1.0 |
| bonobos | 6 | 49.3 | 0.341 | 0.228 ± 0.045 | +2.5 |
| people | 6 | 31.6 | 0.664 | 0.290 ± 0.074 | +5.1 |

The mouse hierarchy shows no excess cyclicity over a matched Bradley–Terry model — if anything it
sits 1.6σ *below* the floor, which deserves a word: a sub-floor value usually indicates the null is
over-dispersed, most often because the fitted rank spacing is too compressed, so the floor itself
may be mis-calibrated downward. Read the mouse result as "no detectable excess" rather than as
positive evidence of unusual transitivity. Subject to that, it sits slightly below the floor — so its antisymmetric social coupling is consistent with a pure
gradient, with no circulating component detectable. That points the same way as the colloid result,
by a route that does not depend on it — **antisymmetry does not imply circulation** — and in a
system where the decomposition is exact rather than inferred. How much weight it can bear is
limited, and the next paragraph is about that limit.

We are deliberately restrained about how much this carries. The two sparse matrices have floors
near 0.58 and cannot resolve cyclicity *in either direction*, so the mouse result should be read as
"consistent with transitive, with low power" rather than as a positive finding — and, as noted
above, a value sitting below its own floor is as likely to indicate a mis-calibrated null as an
unusually transitive group. The two dense
matrices resolve it and disagree with each other, which establishes that the statistic has real
range when the counts support it but says nothing systematic about species. All four are single
unreplicated groups. The binding constraint is interaction density per ordered pair, not data
availability, and the density that round-robin testing protocols produce is an order of magnitude
below what the measurement needs; continuous observation of freely interacting animals is the only
route we know of to close that gap.

A related observation supports the general strategy of this paper rather than any specific claim.
In the longitudinal behavioural table of ref. [27], the *asymmetry* between a directed behaviour and
its reverse — chasing minus being chased, following minus being followed — is a substantially more
reproducible individual trait across days than either direction alone (intraclass correlation
0.576, 0.619, 0.719, 0.781 across four behaviour pairs, against 0.175–0.676 for the components).
The per-group asymmetries also sum to zero to machine precision, which validates the rate columns
internally. The antisymmetric combination being the better-behaved coordinate is exactly the
premise on which the χ construction of Section 2.2 is built.

## 4. Discussion

### 4.1 Isolating a sector requires a parameter that varies it alone

The central methodological result of this work is that testing a claim about the antisymmetric
sector of a coupling requires a parameter that varies that sector and nothing else, and that
obtaining one is harder than it appears. Neither parameter previously proposed for the role
qualifies. The coupling strength α̂ multiplies the symmetric and antisymmetric parts equally and
cancels exactly from their ratio, so it tunes overall interaction strength; sweeping it produced a
trend opposite to the one predicted, for reasons unrelated to reciprocity. The steric size ratio,
as parameterised in the source model, varies packing while the electrohydrodynamic radii remain
pinned. We note that this second failure is specific to the simulation: in experiment a particle's
electrohydrodynamic and steric radii are the same physical size, so the experimental size ratio
does tune nonreciprocity even though the simulated one does not.

Once a genuine knob exists, the causal claim becomes straightforward to establish and turns out to
be large. Nonreciprocity drives the finite-time fragmented morphology called arrested coarsening
(Section 3.2), with a same-particle reciprocal control that the source experiment cannot itself
provide — while, as Sections 3.5 and 4.5 show, simultaneously unjamming the reciprocal gel.

### 4.2 The antisymmetric sector restructures rather than merely circulating

The physical result that most changes the interpretation is that the antisymmetric sector alters
equal-time structure, and does so dominantly. A description that assigns the symmetric part to
energy-like structure selection and the antisymmetric part to circulation cannot be maintained:
both sectors determine the structure, and here the antisymmetric one determines most of it.

This had to be so. The premise that a solenoidal coupling preserves the stationary density holds
only when ∇·(ρ_eq **v**) = 0, a condition that generic nonreciprocal couplings do not satisfy. It
also stands in direct tension with the source experiment, whose central finding — that
nonreciprocity arrests coarsening — is itself a report that nonreciprocity changes steady-state
structure.

Section 3.12 offers a suggestive parallel, which we are careful to label as an analogy rather than
a second instance of the same proposition. In a social dominance network the antisymmetric edge
flow admits an exact, finite-dimensional Hodge decomposition, and the measured flow is consistent
with a pure gradient. But a Hodge curl of a static antisymmetric edge function on a complete graph
of *n* individuals and a probability current in a 2N-dimensional configuration space are different
objects, and Section 1 lists precisely this kind of identification as one to avoid. The dominance
result is therefore a structurally analogous finding in a setting where the decomposition is exact,
not independent confirmation of the colloid result. Given that it is also underpowered
(Section 3.12), we rest nothing on it.

The same decomposition can, however, be computed on the colloid's own contact graph, which is the
directed-graph test the motivating proposal calls for. On every contact edge the antisymmetric
coefficient w_ij = (a_j − a_i)/2 is an antisymmetric edge function: the pair propulsion has
magnitude χ|w| r and points toward the particle with the larger electrohydrodynamic radius.
Because that radius takes two values, w vanishes on L–L and S–S contacts and has one sign on every
L–S contact, so the flow is a two-level potential up to the r-dependence of the kernel. Its
discrete Hodge decomposition on the contact graph (`hodge_contact.py`, six late snapshots per
run) gives a cyclic fraction ‖w − ∇s‖/‖w‖ of 0.077, 0.079, 0.075, 0.071 and 0.063 at χ = 0, 0.25,
1, 1.5 and 8, against 0.73–0.75 for the same magnitudes permuted over edges with random signs:
**the antisymmetric flow is predominantly gradient — approximately 99 % or more in squared norm
at every χ — with a small non-gradient residual** of the size the 11–13 % spread of the kernel
over the contact shell predicts. As an edge function, then, the antisymmetric sector of this
model is close to a two-level potential, with a small residual that is not zero. The caveat of the previous paragraph applies with full force — this is a Hodge
decomposition of a static edge function, not of the configuration-space current — but it is the
colloid-side statement that Section 3.12 could only make by analogy.

Two further reasons no circulation appears deserve to be stated, because they are of different
kinds. The first is symmetry. The force law, the noise and the geometry are all invariant under
reflection, so every spatial pseudoscalar — the angular velocity of a cluster, the vorticity of
the coarse velocity field — has zero ensemble mean at every χ by parity, not by measurement. The
observable-plane signed-area rates of Section 2.4 are not pseudoscalars under reflection and are
not covered by this argument; they remain empirical nulls, and the chiral Vicsek control of
Section 3.11 circulates precisely because it breaks the parity that this model keeps. The second
is empirical and stronger than the four-plane test of Section 3.4. On the dense continuations,
the antisymmetric part of the lagged cross-correlation matrix of seven scalar observables
(cluster count, largest-cluster fraction, mean fragment size, contact-edge count, potential
energy, configurational excess dissipation, edge turnover), summed over lags of 50 to 2000
against a block-shuffled surrogate, gives z = −0.3 ± 1.7, −0.2 ± 1.5, −0.6 ± 0.5, −0.6 ± 1.2,
+0.2 ± 1.1 and −0.8 ± 1.6 at χ = 0, 0.25, 0.5, 0.75, 1 and 1.5 (three runs each; largest single
|z| = 2.2 of eighteen). No pair of these observables traces a cycle at any lag. And in the
cluster-size coordinate itself, what the condensate sheds and what it reabsorbs have the same
size distribution at all six χ (two-sample KS p ≥ 0.53). That is a marginal-distribution result:
it does not test the full transition-current balance of the loop, and it cannot establish that the
loop carries no circulating current, only that the one-dimensional projection onto fragment size
shows none. "No circulation in any projection examined" therefore stands with three more
projections examined, one of them the natural coordinate of the structural transition.

What the antisymmetric sector does uniquely contribute is a positive excess nonconservative work
rate above a structural crossover, quadratic in χ there and unresolved below (Section 3.10). The
configurational decomposition cancels 87–91 % of the unopposed contribution at every drive above
the crossover and, within error, all of it below, so the sector's dissipative content is the
residue the symmetric forces do not hold, and that
residue is naturally represented by a dissipation functional. We put it no more strongly than
that: path actions and dissipation functions are not mutually exclusive descriptions, and the
measurement does not exclude an action representation. What it does exclude is the specific claim
under test — that the antisymmetric sector merely adds a structure-preserving term to an
equilibrium action. We attach no stronger claim than this: Section 3.10 also reports an apparent
antisymmetric cross-coupling to a second nonreciprocal channel, which a structurally dissimilar
control shows not to survive, and we do not rest any conclusion on it.

### 4.3 Network observables are largely blind to nonreciprocity

Several standard network measures fail to detect nonreciprocity here, and they fail in
instructive ways.

Newman modularity is degenerate on these contact graphs. Cluster count varies by a factor of forty
across the pilot-scale conditions while modularity varies by 15%, and the two are uncorrelated; a
three-cluster and a thirty-three-cluster configuration are indistinguishable. The degeneracy
survives a fourfold increase in system size, so it is a property of modularity on two-dimensional
contact networks at this density rather than a finite-size effect. Modularity excess over a
degree-preserving null moreover *decreases* with nonreciprocity and is highest in the reciprocal
reference, inverting the predicted ordering.

Edge turnover fails for four independent reasons: about 70% of its value at χ = 1 survives at χ = 0,
where the dynamics obey detailed balance and the stationary current vanishes; it correlates strongly with cluster morphology; its
lag-one value is dominated by flicker at the contact threshold rather than by rewiring; and it
becomes non-monotonic in χ at paper scale, a trend reversal that only appeared after
equilibration.

The circulation measure is the most interesting failure, because it turns out to be the right
instrument pointed at the wrong question. It registers nothing in the colloid suspension at either
system size, nor in an all-pairs lag test on seven observables (Section 4.2), yet the same
estimator applied to cilia gives a median |z| of 3.2 (3.0 under an amplitude-preserving iAAFT
null). Circulation is
therefore a poor probe of nonreciprocity. It is not, however, a probe of biological organisation
either: a Vicsek-type ensemble with an imposed single-particle turning rate registers more strongly
than cilia do (Section 3.11). What it detects is the presence of a cyclic mode in the chosen
projection, whether that mode is emergent or imposed.

### 4.4 Two kinds of arrest, distinguished by dissipation

Section 3.5 establishes that the reciprocal reference state is a kinetically arrested gel: pair
wells are 4.8 to 7.0 k_BT deep so bonds do not break thermally, cluster diffusion scales as 1/n so
a thousand-particle cluster moves less than a particle diameter over an entire run, and — decisively
— a condensed configuration evolved at χ = 0 with no activity whatever remains condensed. The
condensate is equilibrium-accessible; the reciprocal dynamics simply cannot reach it.

This carries a caution for the wider literature on this system. Both the reciprocal control and the
monodisperse reference are *arrested*, and they are statistically indistinguishable from one
another. "Arrested coarsening" therefore does not by itself indicate an actively maintained state.
What distinguishes the two regimes is dissipation. At the reciprocal limit the nonconservative
work vanishes by construction (zero total entropy production additionally requires the equilibrium
stationary ensemble, and an aging reciprocal gel can still relax dissipatively); the excess is
unresolved at the lowest tested drive and becomes clearly resolved, and quadratic, across the
crossover. Where an analogy to homeostatic or biological organisation is intended, the transferable
quantity is an irreversibility measure rather than a structural one, since a frozen aggregate and a
dynamically maintained one may be structurally similar while differing absolutely in dissipation.

### 4.5 Arrest without a steady state

The finite-cluster state of this model is real, reproducible and long-lived — it is what one
observes on the timescale of the source experiment — but it is not asymptotic. Continued to twice
the source duration, the same system phase-separates. A variational account that attributes the
finite-cluster state to a steady-state balance between energy and connection cost is therefore
describing a transient rather than an attractor.

The transient is nonetheless long, and its lifetime is measured rather than extrapolated:
the condensation time scales as t_x ∝ N^{1.0} over a fourfold range of N, with the lcf(t) curves
collapsing in t/N, and the mass-weighted growth law S_w ∝ t^{0.5–1.0} is consistent with that
scaling (Section 3.8). Under t_x ∝ N, a suspension of 10⁸ particles at the source density would
stay in the finite-cluster regime of order 10⁴ times longer than the source's 22 000-particle
domain, which is long enough that the distinction between long-lived transient and steady state may be
difficult to draw operationally, and it reconciles our result with the source report, which is
accurate for its own system size and duration.

A related inversion deserves explicit statement, because the phrase "arrested coarsening" invites
the opposite reading. The reciprocal case has by far the *slowest* coarsening (mass-weighted
cluster size growing as t^{0.11}, Section 3.8), and every nonreciprocal case up to χ = 3 grows
with an exponent between 0.4 and 1.0 on the same measure. Nonreciprocity accelerates growth
of the majority phase — it unjams the gel — while sustaining a population of small fragments.
"Arrest" in this model therefore denotes a persistent fragment population, not a slowed majority
phase, and the genuinely arrested state is the reciprocal one.

### 4.6 What the biological comparison does and does not establish

The colloid suspension is nonreciprocal, dissipative and strongly self-organising, and it is
uncontroversially not alive. That alone makes it a useful control, and the comparison below does
not depend on the tier question of Section 3.10 — which, as reported there, our measurements
resolve only as a short-lag agreement window. What the comparison requires is only that the colloid be driven, far from
equilibrium, and demonstrably irreversible at the microscopic level, all of which is established
independently in Section 3.4.

Cilia differ from the colloid on a measurable axis — their irreversibility survives coarse-graining
— and it is tempting to read that as a signature of biological organisation. A synthetic control
shows that it is not. A chiral Vicsek-type ensemble, given a turning rate by construction,
circulates more strongly than cilia do, while a nonreciprocal but non-chiral variant circulates not
at all despite far higher polar order. The discriminating property is therefore the presence of a
cyclic mode in the projection, not biology and not nonreciprocity.

We regard this as the more useful outcome. It converts a claim we could not have defended into a
narrower one that the data support: microscopic irreversibility need not project onto macroscopic
observables, and whether it does is a structural property of the system's collective modes rather
than of its provenance. It also identifies what a genuine biological signature would have to do —
distinguish an emergent limit cycle, such as the coordinated ciliary beat, from an imposed one,
such as a single-particle turning rate — which the area-rate estimator by construction cannot.

### 4.7 Consequences for a network-weighted action principle, and for minimum-action learning

One specific formulation of the decomposition tested above is the Network-Weighted Action Principle
(NWAP) [9, 33], which posits that a system minimises

    S_NW = ∫ (E − I + A·C) dt

with E an energetic or operational cost, I the system's information or diversity, A a connectivity
weight and C the cost of forming and maintaining connections. It has been applied to physiological
causation [9], symbolic law discovery [31], neural architecture [32], and marine metabolic networks
[30]. What this study tests is a *directed-graph extension* of it — unpublished, and formulated by
the present author while designing this study — that makes the symmetric/antisymmetric assignment
explicit, together with one published prediction, the modularity-excess signature of ref. [30].

*That framework is the present author's own, and this study was undertaken to test its predictions;
we state the outcome directly and note the conflict of interest.* The outcome is mixed, and the
parts that fail and the parts that survive are informative in different ways.

**The decomposition proved experimentally useful.** The directed-graph extension holds that a
nonreciprocal coupling should be split into symmetric and antisymmetric sectors, with the
antisymmetric sector carrying the nonequilibrium content. Constructing that split is what made
every subsequent measurement possible, and the prediction that cluster statistics are controlled by
the antisymmetric-to-symmetric ratio is borne out across two system sizes, with a same-particle
reciprocal control the source experiment cannot itself provide. We are deliberate about what this
does and does not establish. The decomposition is an algebraic operation, and showing that cluster
statistics depend strongly on χ demonstrates that the antisymmetric sector matters. It does not
validate the NWAP functional, its information term, its connection-cost term, or any predicted
functional relationship between them. The useful thing that survives is the instrument, not the
theory behind it.

We state this more carefully than a rank correlation would suggest. The Spearman coefficient across
the χ sweep is ρ = +0.931, but the underlying relation is *not monotonic*: cluster count dips to 2.2
at χ = 0.25, a factor of five below its reciprocal value, before rising by two orders of magnitude
(Section 3.5). A rank correlation is a poor summary of such a curve. The defensible statement is
that nonreciprocity strongly controls cluster statistics, that the control is non-monotonic, and
that the form which does predict the dip is not a functional at all but a balance of two kinetic
rates with different thresholds (Section 3.5) — which no proposed functional form, NWAP's
included, contains. This is, to our knowledge, the first
NWAP prediction tested in a controlled physical system rather than by meta-analysis, and on this
point it succeeds.

**The mechanism attributed to that sector does not survive.** Three specific claims fail:

*Solenoidality.* The premise that the antisymmetric sector cannot alter a static structural
observable is contradicted directly: cluster count, a single-snapshot quantity, changes by a factor
of 12.4 when the antisymmetric coupling alone is scaled. The premise holds only where
∇·(ρ_eq **v**) = 0, which generic nonreciprocal couplings do not satisfy. It is also in tension
with the source experiment, whose central result — nonreciprocity arrests coarsening — is itself a
statement that nonreciprocity changes steady-state structure.

*Circulation.* No circulation is detectable in any coarse network observable, at either system
size, against a reciprocal null whose stationary current vanishes. The "circulating partition" half of the prediction has
no empirical support in this system, notwithstanding that the system is genuinely irreversible.

*The modularity signature.* This is the one prediction with a published pedigree. Ref. [30]
establishes, in marine metabolic networks, that modularity *excess* over a null — rather than
absolute modularity, which there arises from sparsity alone — is the informative signature of
constrained organisation, reporting ΔQ ≈ 0.15–0.40 over configuration-model, label-permutation and
bipartite-incidence nulls. Carried over to the present system, the directed-graph extension
predicts a sustained excess in the nonreciprocal case relative to a reciprocal null. **It fails
with inverted sign.** Our excess over a degree-preserving (configuration-model) null is of
comparable magnitude to theirs, 0.38–0.46, so the quantity is not simply absent; but it *decreases*
monotonically with nonreciprocity and is largest in the reciprocal reference. The underlying reason
is that Newman modularity is degenerate on these contact graphs — cluster count varies forty-fold
while Q varies by 15%, and the two are uncorrelated — so the test could not have discriminated in
either direction. We note that this diagnosis is specific to two-dimensional contact networks at
this density and does not transfer to the sparse bipartite metabolic networks of ref. [30], where
the excess is measured against a different null on a different topology.

**What the antisymmetric sector actually contributes is dissipation, above a threshold.**
Section 3.10 establishes that it produces a positive entropy production, quadratic in the coupling
once the contact network yields and unresolved before, with 87–91 % of the unopposed
contribution cancelled by the symmetric forces throughout. A positive quadratic in the drive is the functional form of a
dissipation functional, not of an action extremum. An antisymmetric, hence non-dissipative, cross-coupling to a second nonreciprocal channel does
not survive its structurally dissimilar control (Section 3.10), so the claim here rests on the
dissipation alone.

**Where this leaves the framework.** We would draw three conclusions, offered constructively.

First, the broadening that NWAP requires is not ad hoc. There is a standard hierarchy of
variational principles for stochastic dynamics, and least action does not in fact fail for
nonreciprocal systems: the Onsager–Machlup action is well defined for any drift field. What fails
is the equilibrium corollary. NWAP's Triple-Action, with its explicit information term traded
against an energy term under structural constraints, has the formal shape of a Maximum Caliber
functional [10] — the path-entropy principle that reduces to MaxEnt statically and to Onsager near
equilibrium. Placing NWAP as a constrained member of that family, with network-structural rather
than thermodynamic constraints, would give it a principled home and two checkable reductions: it
should reduce to a Rayleighian at weak drive and to free-energy minimisation at zero drive.

This also sharpens NWAP's relation to its neighbours. It has been positioned alongside the
free-energy principle [34] and dissipative adaptation [35], both of which are likewise
functional-minimisation accounts of organisation. The distinction our results suggest is not which
functional is minimised but *what the functional is allowed to govern*: on the evidence here, a
variational account of this kind can describe the dynamics at fixed structure, and the transition
that selects the structure is a separate problem. That boundary applies to the neighbouring
frameworks as much as to NWAP, and is a more useful place to look for discriminating tests than the
choice of functional.

Second, the discipline that keeps such a broadening falsifiable is that each tier carries its own
experimental signature and membership must be *measured*. We tested three such signatures here.
A framework asserting that "some variational principle applies" forbids nothing; one asserting
tier-II membership predicts quadratic dissipation, Onsager reciprocity, and a single effective
temperature, each of which can fail independently. In this system the first is near-tautological
at fixed structure, the second does not survive a structurally dissimilar control, and the third
passes over a short-lag window and fails beyond it. That is a harder outcome to accommodate than a
clean pass or a clean fail, and it is one that no amount of
reinterpretation could have produced from the framework alone. It also suggests that the tier
taxonomy, useful as a way of organising what to measure, is not straightforwardly decidable by
measurement in a system like this one.

Third, the boundary of applicability should be stated rather than obscured, and it can now be
stated in two measured forms. In time: fluctuation–dissipation with a single shared temperature
holds over lags shorter than t_FDT(χ) and fails beyond, with t_FDT falling roughly as χ⁻²
(Section 3.10), so a near-equilibrium description covers contact- and intra-cluster-scale
dynamics and not cluster-scale motion. In structure: the selected state is an approximately stationary kinetic balance in which each
channel pair is approximately balanced, with the non-monotonic fragment count accompanied by
distinct drive dependences of condensate shedding and reabsorption (Section 3.5), and the crossover that
separates fragments from condensate has a horizon-converged committor of one half (Section 3.9). Whatever variational description
turns out to apply — and Section 3.10 does not establish that a Rayleighian does — its scope can
extend at most to the dynamics at fixed structure within that lag window. Structure selection
is not explained by the equilibrium-like or Rayleighian constructions tested here. We put it that way rather than claiming no variational principle covers
it: nonequilibrium quasipotentials, large-deviation theory, macroscopic fluctuation theory and
model-specific constructions all address aspects of structure selection, and we have not tested
them. Since
structure selection is precisely what NWAP was built to explain, this identifies the gap that new
theory would have to fill, and it is a more useful result for the framework than a claim of
coverage would have been.

**Connection to minimum-action learning.** The same programme has a computational arm,
Minimum-Action Learning [31], whose discriminating criterion is energy conservation under dynamical
rollout. That criterion presupposes a conserved energy. Overdamped dynamics have none: at χ = 0 a
potential U exists but is not conserved, since even without noise dU/dt = −Σ_i μ_i|∇_iU|² ≤ 0, and
with noise it fluctuates and exchanges with the bath. What χ = 0 does provide is a potential
whose dissipation is entirely the relaxation −dU/dt, whereas above χ = 0 there is a second,
nonconservative contribution — the excess of Section 3.4 — that no potential accounts for. A
criterion built on a conserved energy therefore cannot be applied as stated; one built on
consistency with the measured excess dissipation could be, and the trajectories generated here
are a ready benchmark for it. We test none of this and note it only as an implication.

**A note on the empirical hook.** Computing modularity excess on published trajectory data,
the reanalysis that would most directly connect NWAP to this literature, fails with inverted
sign for reasons intrinsic to the measure (Section 3.3), and we do not recommend it. The analysis that
does work, and that we would suggest in its place, is the χ decomposition with cluster statistics
and entropy production as the observables, together with the two-kinds-of-arrest distinction of
Section 4.4 as the bridge to the biological claim.

## 5. Limitations

**The reciprocity parameter is synthetic.** χ is a control parameter of our construction rather
than an experimental one, and χ > 1 lies outside the derived electrohydrodynamic model. Values
above unity are reported as trend evidence only. The corresponding physical knob is the
electrohydrodynamic radius contrast, which is experimentally accessible but which we did not sweep
because it confounds the nonreciprocity ratio with the overall coupling strength of the small–small
pair interaction.

**The excess dissipation is a residual of two cancelling terms.** The configurational estimator
avoids the trajectory estimator's discretisation problem, but the excess it returns is the ten
per cent of $\tilde a$ that contact screening leaves, so a shift of a few per cent in $\tilde b$ moves it by a
large fraction at low drive. The production-timestep configurations carry such a shift, and every low-drive number in this
paper is taken from configurations re-equilibrated at dt = 0.005, where the stationarity
residual is 8 % of its production-timestep value. Residual bias of order that 8 % of the correction
cannot be excluded, and the χ = 0.5 entry in particular should be read as "a third of the plateau,
give or take a third of itself"; at dt = 0.0025 on one seed the χ ≤ 0.75 values agree with the
dt = 0.005 values within their statistical errors (Section 3.4), which are now the larger
uncertainty. Below the crossover the excess is unresolved, not shown to vanish. The trajectory probe, where we use it, remains a subtracted
artifact with the caveats stated in Section 2.4.

**Tier membership is established only as a lag-window boundary.** Of the three near-equilibrium
diagnostics tested in Section 3.10, one is near-tautological at fixed structure, one does not
survive a structurally-dissimilar control, and the third — fluctuation–dissipation on
single-particle coordinates — holds with a single shared temperature over lags shorter than
t_FDT(χ) and fails beyond, as a time-domain statement that has not been resolved in frequency. We make no claim about which variational tier the system as a whole occupies; the claim is
that the near-equilibrium description has a measured lag window of validity that shrinks with
drive.
The Harada–Sasa sum rule, which would tie the frequency-resolved violation to the dissipation of
Section 3.4 quantitatively, was not evaluated, and the Green–Kubo prediction of the low-drive
curvature from the χ = 0 reference was not computed.

**Response coefficients require a plateau test.** Both the Onsager and effective-temperature
measurements produce stable, small-error-bar values at any single measurement window or
perturbation strength, and only a scan against the window or the strength shows whether the value
is a coefficient or an artefact of the protocol. We therefore distinguish response estimates that have reached a plateau from measurements
that remain window-dependent — the dissimilar-channel cross-response of Section 3.10 is one —
and do not interpret the latter as asymptotic coefficients. Readers should treat any such
coefficient reported without that check, in this literature generally, with corresponding
caution.

**The critical-size crossover rests on two criteria that do not coincide.** The committor of the
full fragment dynamics is converged in the horizon and crosses one half at $S_0$ ≈ 8, but the
mean-drift criterion of classical nucleation theory is met only within the monomer-exchange
channel (zero at S ≈ 10), and the mean drift of the full process has no zero. Both are measured
at one system size and density, with six parent configurations at the longer horizon.

**The directed contact graph is a static object.** Section 4.2 builds it and finds its
antisymmetric flow approximately 99 % gradient in squared norm, but that is a Hodge decomposition of an edge function at one
instant, not of the configuration-space probability current, and the two are different objects;
the result constrains what the pair law can do, not what the dynamics do.

**Box geometry in the source-specification replicate.** The source domain is 648 × 360 µm, an
aspect ratio of 1.8, whereas the replicate uses an equal-area square. We tested this directly with
a rectangular-box integrator (`simulate_rect.py`, validated as bitwise identical to the square
kernel in the L_x = L_y limit): at N = 4000, χ = 1.5, matched area, density and duration, the 1.8:1
rectangle gives lcf = 0.804 ± 0.008 against 0.779 ± 0.004 for the square, a difference of 3%.
The square approximation is therefore justified at this system size, though we have not repeated
the test at N = 22,000.

**σ²_v is not comparable to the source definition.** The source reports 5 × 10⁻³ and 2 × 10⁻²
µm²/s² for monodisperse and bidisperse respectively, a ratio of 4; ours is approximately 90. Ours
is the variance of *speed* coarse-grained over one snapshot interval rather than an instantaneous
velocity variance, so the discrepancy is most likely definitional, but this cannot be confirmed
without the source sampling interval.

**Pairwise forces only.** Many-body hydrodynamics are absent, inherited from the source model,
which shifts absolute cluster sizes and may bear directly on the transience result of Section 3.8:
hydrodynamic interactions are exactly the class of ingredient that could stabilise finite clusters
indefinitely. The nucleus of Section 3.9 and the rate balance of Section 3.5 are likewise
properties of the pairwise model.

**The χ = 0 entropy production is zero by construction, not by measurement.** The paired protocol
subtracts each configuration's own χ = 0 bias probe, so at χ = 0 the estimate and its control are
the same quantity and the net is identically zero. The χ = 0 row therefore validates nothing; the
argument that χ = 0 has zero current is the analytic one of Section 2.3, not an empirical one. What
the paired protocol does establish is that the *excess* over that reference is resolved for
χ ≳ 0.75.

**The dominance analysis is underpowered where it matters most.** Of the four published matrices in
Section 3.12, the two built by round-robin testing have 3.5 and 4.3 interactions per ordered pair
against a matched transitive floor near 0.58, so they cannot resolve cyclicity in either direction.
The rodent result is therefore "consistent with a gradient, with low power" and not a positive
finding. All four are single unreplicated groups, and the floor is computed under a Bradley–Terry
null; a different transitive generative model would shift it.

**Sample sizes in the biological comparison.** The cilia analysis uses two shape modes. Under the
stricter iterative amplitude-adjusted Fourier-transform surrogate, which preserves the amplitude
distribution as well as the power spectrum, the pooled median |z| falls from 3.2 to 3.0 and the
fraction above |z| = 3 from 0.55 to 0.51, with 92 % above |z| = 2 under both nulls; the result is
not sensitive to the choice. Beat frequency is a spectral-peak estimate and is not monotonic at
the top of the ATP range.

### 5.1 Pre-registration against outcome

Thresholds for the χ experiment were committed in `EXPERIMENT.md` before any χ run was executed.
Setting the pre-registered predictions against what happened: two of the four hypotheses failed,
and one of the failures ruled out the study's original headline interpretation.

| # | pre-registered prediction | pre-committed threshold | outcome |
|---|---|---|---|
| H1 | edge turnover increases monotonically in χ | Spearman ρ > 0, p < 0.01 | **passed at pilot scale, failed at paper scale.** Turnover is non-monotonic in χ, peaking at χ = 0.75 (§3.7). The monotonicity was an artifact of not having reached steady state. |
| H2 | structure moves far less than currents ("frozen Q, circulating partition") | R_Q < 0.1 **and** R_turnover > 0.5 | **failed as intended, and informatively.** Q did stay within threshold, but not because structure is preserved: n_cl moves by 12.4× (§3.2) and Q is degenerate on these graphs (§3.3). The pre-committed reading "turnover rises, Q flat ⇒ prediction confirmed" would have been wrong, which is why §3.3 tests Q's discriminating power directly. |
| H3 | σ²_v increases with χ | directional | **passed.** ρ = +0.975 (nominal p = 7×10⁻¹², treating blocked runs as independent; §3.2). |
| H4 | at χ = 0 turnover falls toward the monodisperse baseline | if not, mono/bi gap was polydispersity | **partial, and the informative arm.** Turnover at χ = 0 is 0.260 against 0.148 monodisperse and 0.378 at χ = 1, so about half the gap is polydispersity (§3.6). The pre-committed reading — that this "falsifies the project's headline interpretation while leaving the measurement intact" — is what happened, and is why edge turnover is retired as a probe in §3.7. |

Two points about this table. First, the one hypothesis that passed cleanly (H3) is the one whose
observable, σ²_v, survived scrutiny; the two that failed are the ones built on edge turnover and
modularity, both of which §3.7 and §3.3 go on to disqualify. Second, the pre-committed failure
readings were followed rather than reinterpreted: H4's stated consequence was that the headline
interpretation dies, and it did. A design whose failure modes were committed in advance and then occurred is the
pre-registration working as intended, not evidence against the design.

## 6. Conclusions

1. Isolating the antisymmetric sector of a nonreciprocal coupling requires a parameter that varies
   it alone. Neither the coupling strength nor the steric size ratio does so; the reciprocity
   mixing parameter χ does. At χ = 0 the *dynamics* are equilibrium-compatible — conservative forces
   with a uniform fluctuation–dissipation ratio, hence detailed balance and zero stationary current
   — which is a stronger reference than a control, though the finite-time configurations it
   realises are kinetically trapped and do not sample the corresponding Boltzmann distribution.

2. Nonreciprocity produces the finite-time fragmented morphology conventionally called arrested
   coarsening, in the operational sense defined in Section 3.2 — many coexisting, continuously
   reorganising clusters — while simultaneously unjamming the reciprocal gel. At N = 4000 the
   cluster count rises by a factor of 12.4 from χ = 0 to χ = 1.5 with the reciprocal coupling held
   bit-for-bit fixed; the pilot scale gives 9.8× for the same contrast, so the effect is reproduced
   at a second system size but the factor itself is the N = 4000 number.

3. The antisymmetric sector alters *static* structure and is not structure-preserving. The
   solenoidal premise fails, and would have had to, since it contradicts the source experiment's own
   central finding.

4. Newman modularity is degenerate on these networks and its excess over a degree-preserving null
   runs opposite to the predicted direction. Cluster-size statistics, not modularity, carry the
   signature.

5. No circulation is detectable in the coarse network observables examined — nor in an all-pairs
   lag test on seven observables, nor in the cluster-size coordinate, nor, by parity, in any
   spatial pseudoscalar — yet the excess dissipation over the reciprocal control is positive and
   resolved above a crossover. It is unresolved for χ ≤ 0.25, a third of the quadratic plateau at
   χ = 0.5 and on it from χ = 0.75 to 8; the configurational decomposition cancels 87–91 % of the
   unopposed antisymmetric contribution above the crossover, through the steric contact force,
   and all of it within error below, through the pair attraction. An apparent negative linear term in a two-term fit is this cancellation, quadratic and
   structural. The dynamics are genuinely irreversible; that irreversibility did not appear in
   any projection we examined, which is a statement about those projections rather than about
   the dynamics, and the antisymmetric pair law is itself predominantly gradient on the contact
   graph, with a small non-gradient residual.

6. The reciprocal reference state is a kinetically arrested gel, not an equilibrium structure. Weak
   nonreciprocity unjams it, producing a non-monotonic dependence of cluster count on χ, and the
   genuinely arrested state is the reciprocal one.

7. No finite characteristic cluster size is selected: the largest cluster grows as N^1.00 ± 0.07
   across a fourfold range of system size, and the largest-cluster fraction is flat (0.79, range
   0.71–0.86) across that range without any fit. The scale that is selected, S* ≈ 8–10 particles,
   is a **critical-size crossover**: the committor of the fragment-size dynamics, converged in the
   observation horizon, crosses one half at $S_0$ = 7.9 (7.6–8.3 over parent configurations; 8.4
   without absorption), and the monomer-exchange drift changes sign at 10.1 (9.7–10.7). On the
   committor evidence it is a critical nucleus of the fragment dynamics; the mean-drift criterion
   holds only within the monomer-exchange channel. It is independent of N, of density within a
   band, and of χ above the crossover. The finite-cluster
   state is a long-lived transient whose lifetime is measured, not extrapolated: the condensation
   time scales as N^{1.0} and the largest-cluster-fraction curves collapse in t/N over a fourfold
   range of N.

8. **The tested single-particle fluctuation–response relation has a measured short-lag
   agreement window, not a tier.** Of the three diagnostics
   we tested, entropy production is quadratic above the structural threshold but that form is
   close to guaranteed by the construction of χ; a differential cross-response between two
   nonreciprocal channels is predominantly antisymmetric when the channels share a structure, but
   a structurally dissimilar control inverts the pattern, so the reciprocity reading does not survive.
   The third, measured on single-particle coordinates with a random-sign force conjugate to the
   single-particle fluctuation, calibrates to unity within 3 % at every lag at χ = 0 and, out of
   equilibrium, gives a single effective temperature equal to T shared by both species over lags shorter
   than t_FDT(χ) — 216, 37 and 13 at χ = 0.5, 1 and 1.5 — with a shared violation growing as t
   at longer lags; this is a time-domain statement, not a frequency-resolved one. On species coordinates the same
   question cannot be posed, because the two coordinates are one degree of freedom up to sign. Whatever variational description
   applies, its scope is at most the dynamics at fixed structure within that lag window, since
   structure selection is a kinetic balance that no functional tested here reproduces.

9. The difficulty in conclusion 8 is itself part of the finding: of three commonly invoked
   near-equilibrium signatures, two are defeated by their own construction rather than by the
   physics — one near-tautological given how χ is built, one an artifact of giving two response
   channels the same structure — and the third gave opposite verdicts on two coordinates until the
   conjugate pair was chosen correctly. Tier taxonomies are useful for deciding what to measure and
   harder than they appear to decide by measurement. Chasing the third on species coordinates also
   turned up a mechanism: nonreciprocal forces do not sum to zero, so the centre of mass drifts,
   near-ballistically over short windows (t^1.96 against t^1.01 at χ = 0). That drift is
   finite-size — amplitude falling roughly as N^−1/2 — and decorrelates over long runs, so we
   report it as a bounded observation and not as a new phenomenon. The colloid suspension none the
   less serves as a control for systems that are driven but not alive, which requires only that it
   be far from equilibrium and demonstrably irreversible — both established independently.

10. For the Network-Weighted Action Principle specifically, the decomposition it proposes
    provided a useful experimental control, and cluster statistics do depend strongly, though
    non-monotonically, on the antisymmetric-to-symmetric ratio; these results do not validate the
    NWAP functional, and the specific mechanism attributed to the antisymmetric sector —
    solenoidality, circulation, and a modularity signature — is unsupported on all three counts. The sector's measurable contribution is a positive excess dissipation above a structural
    crossover, quadratic there and unresolved below, with 87–91 % of the unopposed contribution
    cancelled by the symmetric forces throughout. That form is naturally represented by a dissipation functional; it does not
    support the stronger claim that the antisymmetric sector merely adds a structure-preserving
    term to an equilibrium action. The structure the sector selects is an approximately stationary kinetic balance in which
    shedding and reabsorption, fission and fusion, and nucleation and dissolution are each
    approximately balanced, which is the concrete thing a functional for structure selection
    would have to reproduce.

11. In published social-dominance matrices, where the symmetric/antisymmetric decomposition is
    exact and finite-dimensional rather than inferred, the two sparsely sampled matrices (mouse,
    caribou) are consistent with a pure gradient at an interaction density where the measurement
    has limited power, while the two densely sampled ones (bonobos, people) show a cyclic excess
    over their transitive floor at 2.5σ and 5.1σ. The dominance evidence is therefore mixed, and
    "consistent with a gradient" applies to the underpowered pair only; it is not an independent
    confirmation of the colloid result.
    On the colloid's own contact graph the same decomposition is exact and gives a predominantly
    gradient field, approximately 99 % in squared norm, at every χ. Separately, the Hodge decomposition of dominance
    must be computed on log-odds: on raw counts a perfectly transitive group reads as 16%
    intransitive however much data is collected, which is a floor that published intransitivity
    claims should be checked against.

12. Cilia show macroscopic circulation in shape space where the colloid, at both system-averaged
    and single-cluster level, shows none. A synthetic chiral flock, however, circulates more
    strongly than cilia while a nonreciprocal non-chiral flock shows none at all, so the measured
    property is the presence of a cyclic collective mode rather than biological organisation. The
    defensible conclusion is that microscopic irreversibility need not project onto macroscopic
    observables, and that whether it does is a structural property of the system's collective
    modes.

## Figures

Generated by `make_figures.py` from the analysis outputs; no values are hand-entered.

**Figure 1 — The χ construction and its validation.** (a) Force coefficients c_i (red) and c_j
(blue) against pair separation at χ = 0, 0.5, 1. At χ = 0 the two coincide, restoring Newton's
third law exactly. (b) The symmetric sector (blue, five superposed curves at
χ = 0, 0.25, 0.5, 1, 1.5) is invariant in χ, while the antisymmetric sector (red) is exactly
linear in it. (c) The antisymmetric-to-symmetric ratio for four values of the coupling strength
α̂; the curves are indistinguishable, showing that α̂ cancels identically from the ratio and
therefore cannot serve as a reciprocity knob. `fig1_construction.png`

**Figure 2 — Structure against χ at paper scale (N = 4000, three seeds).** (a) Cluster count, log
axis, showing the non-monotonic dip at χ = 0.25 to a value five times *below* the reciprocal case
before rising by two orders of magnitude; dashed line is the monodisperse reference. (b)
Largest-cluster fraction. (c) Activity σ²_v (orange, log axis) against Newman modularity (blue,
right axis) on a scale spanning only 0.86–0.95: the currents rise by two orders of magnitude while
Q barely moves. Error bars are the standard error over seeds. `fig2_structure.png`

**Figure 3 — Box scaling at fixed density.** (a) Largest cluster against N for three converged
system sizes, with the phase-separation (N¹) and finite-characteristic-size (N⁰) expectations. The
measured exponent is 1.00 ± 0.07. (b) The largest-cluster fraction is flat across a fourfold range of N,
which is the same statement without a fit. `fig3_scaling.png`

**Figure 4 — The critical-size crossover.** (a) Channel-resolved size drift ⟨ΔS⟩ per interval
against cluster size at χ = 1.5, from monomer and dimer exchange events only (Δt = 50, three
seeds, 209,748 cluster observations), with the sign change at S* = 10.1 (68 % bootstrap
9.7–10.7). (b) Committor q($S_0$) of the full fragment-size dynamics from 72 relaunches of
t̂ = 3000 with fresh noise (5,664 fragment fates, 93 % decided near the crossover), with and
without absorption into the condensate counted as the condensed outcome, crossing one half at
$S_0$ = 7.9 (8.4); the dotted curve is the fraction still undecided at t̂ = 3000. `fig4_rates.png`

**Figure 5 — Dissipation and response.** (a) Excess dissipation per particle divided by χ²,
against χ (ordinate in units of 10⁻³), from the configurational estimator at N = 4000; the shaded
band is the plateau of 2.3–3.3 × 10⁻³ and the shaded interval in χ the crossover of Section 3.4.
(b) Symmetric and antisymmetric components of the differential cross-response against
measurement-window duration, for structurally identical (filled) and dissimilar (open)
perturbation channels. (c) The single-particle fluctuation–response ratio, expressed as an
effective temperature relative to the bath temperature, against lag, for the large (solid) and
small (dashed) species at χ = 0, 0.5, 1 and 1.5. Agreement with the equilibrium value over short
lags is an operational time-domain result, not a frequency-resolved test of
fluctuation–dissipation relations. `fig5_thermo.png`

**Figure 6 — Circulation at matched level of description.** (a) Median |z| of the signed area rate
for a single tracked colloidal cluster at three values of χ, against a single *Chlamydomonas*
axoneme, using the identical estimator, shape-descriptor reduction to two modes, and
phase-randomised surrogate protocol. (b) The fraction of objects exceeding |z| = 2. Single active
colloidal clusters show no circulation; axonemes show it in 92% of cases. `fig6_circulation.png`

Supplementary figures: `modularity_test.png` (baseline four-panel), `snapshots.png` (final-frame
renders), `coarsening.png` (trajectories by χ), `chi_test.png`, and `sweep.png` (the
coupling-strength and size-ratio sweeps of Section 3.1).

## Data and code availability

All code, per-run summaries and figures are in the repository
`github.com/martinfrasch/nonrecip-modularity`. Simulation: `simulate.py` (model,
reciprocity parameter, sweeps), `simulate_par.py` (thread-parallel force kernel, validated to
8.9 × 10⁻¹⁵ against the serial kernel), `simulate_rect.py` (rectangular box, bitwise identical to
the square kernel at L_x = L_y), `continue_run.py` (run continuation). Analysis: `analyze.py`
(network observables), `epr_probe.py` (entropy production), `onsager.py`, `onsager_equil.py` and
`onsager_equil2.py` (cross-response), `fdt.py` (effective temperature), `arrest_test.py`
(configuration-swap test), `crossover.py` (per-cluster propulsive force), `irreversibility.py` and
`benchmark_irrev.py` (time-series estimators and their calibration against an exactly solvable
biased ring walk). Comparison systems: `cilia_analysis.py`, `single_cluster_circulation.py`,
`vicsek.py`, `hodge_dominance.py`, `forkosh_asymmetry.py`. Controls:
`drift_control.py` (whether the collective-coordinate superdiffusion is whole-system drift) and
`fdt_driftfree.py` (effective temperature with and without that drift, from the same starts).
Kinetics, fates, response and dissipation: `static_epr.py` (configurational excess dissipation with
small-timestep calibration), `cluster_kinetics.py` (event classification, Becker–Döring drift,
condensate shedding and reabsorption), `committor.py`, `hs_fdt.py` and `hs_analysis.py`
(random-sign fluctuation–response twins and the single-particle effective temperature),
`coarsening_law.py` (structure-factor length, mass-weighted cluster size, chained continuations),
`hodge_contact.py` (directed contact graph and its Hodge decomposition), `lag_asymmetry.py`
(all-pairs lag asymmetry against block-shuffled surrogates); `OPEN_QUESTIONS_PLAN.md` and
`RESULTS_OPEN_QUESTIONS.md` record the programme and its numbers.
Figures: `make_figures.py`. Validation: `tests/test_reciprocity.py` (the χ construction).

Supplementary figures, not reproduced here: `modularity_test.png` (baseline four-panel),
`snapshots.png` (final-frame renders), `coarsening.png` (trajectories by χ), `chi_test.png`, and
`sweep.png` (the coupling-strength and size-ratio sweeps of Section 3.1).

Supporting documents: `AUDIT.md` (observable validity), `EXPERIMENT.md` (pre-registered protocol
with thresholds committed before running), `BOX_SCALING.md`, `COARSENING.md`, `ONSAGER_RESULT.md`,
`CILIA_RESULT.md`, `DOMINANCE_HODGE_RESULT.md`, `FORKOSH_RESULT.md`, `VARIATIONAL_FAMILY.md`,
`FOLLOWUP.md`, `TODO.md`.

Raw trajectories (approximately 4 GB) are not version-controlled. The cilia dataset is openly
available from Dryad under CC0 [6].

## References

[1] S. Hara, M. Okada, K. Kittaka, S. Tanami, Y. Iwasaki, H. Ishikawa, K. Yoshii, Y. Sumino.
*Arrested coarsening in active colloidal suspensions driven by nonreciprocal electrohydrodynamic
interactions.* Phys. Rev. Lett. **137**, 068302 (2026).
DOI 10.1103/96ky-d1p9; arXiv:2509.23164.

[2] M. Fruchart, R. Hanai, P. B. Littlewood, V. Vitelli. *Non-reciprocal phase transitions.*
Nature **592**, 363 (2021).

[3] M. E. J. Newman. *Modularity and community structure in networks.* PNAS **103**, 8577 (2006).

[4] S. Fortunato, M. Barthélemy. *Resolution limit in community detection.* PNAS **104**, 36 (2007).

[5] K. Sekimoto. *Stochastic Energetics.* Lecture Notes in Physics 799, Springer (2010).

[6] V. F. Geyer, J. Howard, P. Sartori. *Ciliary beating patterns map onto a low-dimensional
behavioural space.* Nature Physics **18**, 1465 (2022), DOI 10.1038/s41567-021-01446-2. Dataset:
Dryad, DOI 10.5061/dryad.0gb5mkm2j.

[7] L. Onsager. *Reciprocal relations in irreversible processes.* Phys. Rev. **37**, 405 (1931).

[8] H. B. G. Casimir. *On Onsager's principle of microscopic reversibility.* Rev. Mod. Phys. **17**,
343 (1945).

[9] M. G. Frasch. *Causal thinking in physiology: a search for vertically organising principles.*
The Journal of Physiology (2026), DOI 10.1113/JP290762. *(Defines the Network-Weighted Action
Principle, S_NW = ∫(E − I + A·C) dt. Author's own prior work; see the disclosure in Section 4.7.)*

[10] E. T. Jaynes. *Macroscopic prediction*, in Complex Systems — Operational Approaches, Springer
(1985), for Maximum Caliber; see also P. D. Dixit et al., J. Chem. Phys. **148**, 010901 (2018).

[11] L. F. Cugliandolo. *The effective temperature.* J. Phys. A **44**, 483001 (2011).

[12] C. Battle, C. P. Broedersz, N. Fakhri, V. F. Geyer, J. Howard, C. F. Schmidt,
F. C. MacKintosh. *Broken detailed balance at mesoscopic scales in active biological systems.*
Science **352**, 604 (2016).

[13] F. S. Gnesotto, F. Mura, J. Gladrow, C. P. Broedersz. *Broken detailed balance and
non-equilibrium dynamics in living systems: a review.* Rep. Prog. Phys. **81**, 066601 (2018).

[14] J. Gladrow, C. P. Broedersz, C. F. Schmidt. *Nonequilibrium dynamics of probe filaments in
actin–myosin networks.* Phys. Rev. E **96**, 022408 (2017).

[15] A. V. Ivlev, J. Bartnick, M. Heinen, C.-R. Du, V. Nosenko, H. Löwen. *Statistical mechanics
where Newton's third law is broken.* Phys. Rev. X **5**, 011035 (2015).

[16] S. Saha, J. Agudo-Canalejo, R. Golestanian. *Scalar active mixtures: the non-reciprocal
Cahn–Hilliard model.* Phys. Rev. X **10**, 041009 (2020).

[17] Z. You, A. Baskaran, M. C. Marchetti. *Nonreciprocity as a generic route to traveling states.*
PNAS **117**, 19767 (2020).

[18] S. A. M. Loos, S. H. L. Klapp. *Irreversibility, heat and information flows induced by
non-reciprocal interactions.* New J. Phys. **22**, 123051 (2020).

[19] É. Fodor, C. Nardini, M. E. Cates, J. Tailleur, P. Visco, F. van Wijland. *How far from
equilibrium is active matter?* Phys. Rev. Lett. **117**, 038103 (2016).

[20] C. Nardini, É. Fodor, E. Tjhung, F. van Wijland, J. Tailleur, M. E. Cates. *Entropy production
in field theories without time-reversal symmetry.* Phys. Rev. X **7**, 021007 (2017).

[21] M. Doi. *Onsager's variational principle in soft matter.* J. Phys.: Condens. Matter **23**,
284118 (2011).

[22] P. Sartori, V. F. Geyer, A. Scholich, F. Jülicher, J. Howard. *Dynamic curvature regulation
accounts for the symmetric and asymmetric beats of Chlamydomonas flagella.* eLife **5**, e13258
(2016).

[23] B. H. Good, Y.-A. de Montjoye, A. Clauset. *Performance of modularity maximization in
practical contexts.* Phys. Rev. E **81**, 046106 (2010).

[24] T. P. Peixoto. *Descriptive vs. inferential community detection in networks.* Cambridge
University Press (2023).

[25] R. Kubo. *The fluctuation-dissipation theorem.* Rep. Prog. Phys. **29**, 255 (1966).

[26] U. Seifert. *Stochastic thermodynamics, fluctuation theorems and molecular machines.* Rep.
Prog. Phys. **75**, 126001 (2012).

[27] O. Forkosh, S. Karamihalev, S. Roeh, et al., A. Chen. *Identity domains capture individual
differences from across the behavioral repertoire.* Nature Neuroscience **22**, 2023 (2019).
Data: github.com/OrenForkosh/IdentityDomains (MIT licence).

[28] X. Jiang, L.-H. Lim, Y. Yao, Y. Ye. *Statistical ranking and combinatorial Hodge theory.*
Mathematical Programming **127**, 203 (2011).

[29] J. P. Curley. *compete: analyzing competitive interaction data.* R package (2016),
github.com/jalapic/compete.

[30] M. G. Frasch. *Modularity emerges from action-functional constraints in marine metabolic
networks: a biology-scale validation of the Network-Weighted Action Principle.* arXiv:2605.05254
(2026). *(Author's own prior work: the source of the modularity-excess prediction tested in
Sections 3.3 and 4.7.)*

[31] M. G. Frasch. *Minimum-Action Learning: energy-constrained symbolic model selection for
physical law identification from noisy data.* arXiv:2603.16951 (2026). *(Author's own prior work.)*

[32] M. G. Frasch. *minAction.net: energy-first neural architecture design — from biological
principles to systematic validation.* arXiv:2604.24805 (2026). *(Author's own prior work.)*

[33] M. G. Frasch. *Bridging mathematical formalism and physical law: a network-weighted action
framework for foundation model training.* HAL preprint hal-05347658 (2025). *(Author's own prior
work.)*

[34] K. Friston. *The free-energy principle: a unified brain theory?* Nature Reviews Neuroscience
**11**, 127 (2010).

[35] J. L. England. *Statistical physics of self-replication.* J. Chem. Phys. **139**, 121923
(2013).
