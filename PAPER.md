# Isolating the antisymmetric sector of a nonreciprocal colloidal model: kinetic unjamming, transient arrested coarsening, and hidden irreversibility

**Martin G. Frasch**

*Computational study, August 2026. Repository: `nonrecip-modularity`.*

---

## Abstract

Nonreciprocal interactions do not derive from a scalar potential, which removes the equilibrium
apparatus — Boltzmann statistics, structure selected by energy minimisation, detailed balance —
even though a path-space variational representation survives. A widely held intuition holds that
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

Across 133 simulations at four system sizes we find that the antisymmetric sector alters *static*
structure by an order of magnitude, refuting the solenoidal intuition; that the reciprocal
reference is a kinetically arrested gel which weak nonreciprocity *unjams*, giving a non-monotonic
dependence of cluster count on χ; that the finite-cluster state is a long-lived transient rather
than a steady state, with a crossover time scaling as N^(1/z), z ≈ 0.3; and that no stable
characteristic cluster size is selected. Entropy production is positive and resolved, yet no
circulation is detectable in coarse observables — including at single-cluster level, matched to a
comparison with beating *Chlamydomonas* axonemes, where circulation is strong. Several standard
network measures, Newman modularity among them, prove unable to detect nonreciprocity at all.

**Keywords:** nonreciprocal interactions, active matter, arrested coarsening, entropy production,
kinetic arrest, coarse-graining

---
## 1. Introduction

Nonreciprocal interactions, in which the force that *i* exerts on *j* is not the negative of the
force *j* exerts on *i*, do not derive from a scalar potential. It is worth being precise about
what this does and does not preclude, because the two are often conflated. It does *not* preclude a
variational representation of the dynamics: the Onsager–Machlup action is well defined for any
drift field, gradient or otherwise, and path-space action principles apply to dissipative and
stochastic systems generally. What it precludes is the *equilibrium corollary* — a Boltzmann
stationary distribution, structure selected by minimising an energy, and detailed balance.

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
Variants of this decomposition appear across the nonreciprocal-matter literature [15–18], and one
specific formulation — a network-weighted action principle whose directed-graph extension makes
exactly this assignment [9] — supplied the immediate motivation for this study. We return to that
formulation in Section 4.7; the predictions tested below are consequences of the decomposition
itself and do not depend on it.

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

We test all three predictions. The third is confirmed; the first two are not.

Testing that proposal requires a control parameter that varies the ratio of antisymmetric to
symmetric coupling and nothing else. We show in Section 3.1 that neither parameter previously
proposed for the role does so, and we construct one that does. With that instrument in hand, three
distinct questions become answerable, and they organise the paper.

The first is causal: what does nonreciprocity actually do to this system? Section 3.2 establishes
that it drives arrested coarsening, and that it does so by altering static structure — which
directly contradicts the solenoidal premise. Sections 3.3 to 3.7 then examine the observables
themselves, and find that several standard network measures are unable to detect nonreciprocity
at all. Sections 3.5 and 3.8 to 3.9 characterise the resulting steady state, showing that the
reciprocal reference is a kinetically arrested gel rather than an equilibrium structure, that no
finite characteristic cluster size is selected, and that the scale the system does select is a
nucleation barrier rather than a stable organisational level.

The second question is which variational tier the system occupies. Rather than asserting that some
principle applies, we treat tier membership as an experimental matter with measurable signatures,
and test three of them: quadratic dissipation, Onsager reciprocity, and a single effective
temperature (Section 3.10). The system passes two of the three, while the third proves ill-posed out of equilibrium, which places it in the linear-response tier on every diagnostic that admits a definite answer and identifies a
Rayleighian, not an action, as the appropriate variational object for its dynamics.

The third question follows from the second. If nonreciprocity, dissipation and self-organisation
together are insufficient to leave linear response, they cannot be what distinguishes living
organisation. Section 3.11 therefore uses the colloid suspension as a control — a well
characterised system that is driven but not alive — and compares it against published recordings
of beating cilia. The two differ sharply, and on an axis that neither the modularity measures nor
the entropy production alone would have identified.

We use the source model without modification except where explicitly stated, and we validate our
reimplementation against the published specification and duration in Section 3.8.

## 2. Model and methods

### 2.1 Equations of motion

We reimplemented the model of ref. [1] (End Matter Sec. II, Eqs. 7–9). Overdamped dynamics in
two dimensions with periodic boundaries, in nondimensional units (length λ = 9 μm, time
τ = 0.01 s):

    dx_i/dt = ξ_i(t) + (1/s_i) [ Σ_j F^col_ij + Σ_j F^EHD_ij ]

`F^col` is a reciprocal soft-core linear spring of natural length s_i + s_j. The EHD term is

    F^EHD_ij = −α · l_j⁴ / (r² + l_j²)^{5/2} · **r**_ij ,    r ≤ 1

and is nonreciprocal because the force on *i* scales with *l_j*, the **other** particle's EHD
radius. Noise satisfies ⟨ξ_i ξ_j⟩ = (σ²/s_i) δ_ij δ(t−t′), σ = 2.6×10⁻³. Baseline parameters
follow ref. [1]: α̂ = 0.005, l_L = s_L = 1/6, l_S = s_S = 1/9, steric polydispersity
ε ~ N(0, 1/30), composition N_L : N_S = 1 : 3.

Integration is Euler–Maruyama with dt̂ = 0.05 (stability bound ≈ 0.2), using a hand-rolled cell
list in a single compiled kernel.

### 2.2 The reciprocity mixing parameter χ

Writing g_k ≡ α l_k⁴/(r²+l_k²)^{5/2} and ḡ ≡ (g_i+g_j)/2, we replace the EHD coefficients by

    c_i = (1−χ)·ḡ + χ·g_i ,   c_j = (1−χ)·ḡ + χ·g_j
    F_i = −c_i **r**_ij ,     F_j = +c_j **r**_ij

so that the **symmetric part (c_i+c_j)/2 = ḡ is independent of χ**, while the antisymmetric
part (c_j−c_i)/2 = χ(g_j−g_i)/2 is exactly linear in χ. χ=1 recovers ref. [1]; χ=0 imposes
Newton's third law exactly.

Validation against the compiled kernel (`tests/test_reciprocity.py`): an isolated noise-free
pair shows centre-of-mass drift 0.000×10⁰ at χ=0, 3.95×10⁻² at χ=1, and exactly half that at
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
coupling and not to the reference itself. This is a stronger reference than the monodisperse comparison used
previously, which differs from the bidisperse system in reciprocity *and* polydispersity
simultaneously.

### 2.4 Observables

Contact graph: edge if r < 1.1(s_i+s_j), minimum image. Per snapshot we compute Newman
modularity Q of the Louvain partition, modularity of a degree-preserving (double-edge-swap)
null, the number of connected components with ≥2 particles (n_cl), and the largest-cluster
fraction (lcf). Per consecutive pair we compute the adjusted Rand index between partitions
(reported alongside a same-graph two-seed noise floor), the Jaccard distance between contact
edge sets, and σ²_v.

Two further observables were introduced during this work:

**Circulation.** The signed area rate in a plane of two observables, A = ½⟨x ẏ − y ẋ⟩, is
exactly zero under detailed balance for *any* observable pair and nonzero iff time-reversal
symmetry is broken.

**Entropy production.** The Stratonovich heat Σ_i F_i ∘ dx_i, accumulated by midpoint rule in
the kernel. This vanishes identically under detailed balance and is therefore a current, unlike
turnover or σ²_v. It is unusable at the production timestep — it subtracts two nearly-cancelling
O(μ|F|²) terms and the discretisation residual dominates — so it is measured at dt = 6.25×10⁻⁴
restarting from equilibrated configurations, with a paired same-configuration χ=0 control to
remove each configuration's own bias.

### 2.5 Simulation campaign

129 runs. Pilot scale N=1000, box 12, t̂=1×10⁵ (61 runs: baseline 5 seeds × 2 cases, α̂ sweep
4×3, size-ratio sweep 5×3, χ sweep 6×5). Paper scale N=4000, box 24, t̂=4×10⁵, 200 snapshots
(21 runs: χ sweep 6×3 plus monodisperse reference ×3). Box-scaling series at fixed density
6.944 particles/area: N=9000/L=36 and N=16000/L=48 at χ ∈ {1, 1.5} (6 runs). Replicate at the
source paper's specification: N=22,000, 22.7% type-I, equal-area square L=53.67, t̂=3.6×10⁵
(2 runs). Eleven continuation runs extend selected conditions to t̂=8×10⁵ (or 7.2×10⁵ for the
replicate). Analyses use the late half of each run. Seeds are blocked: seed *k* gives the same
initial configuration and noise stream at every χ.

Tier-membership measurements (Section 3.10) add: entropy production probes at dt = 6.25e-4 from
equilibrated configurations, with a paired same-configuration control; Onsager cross-coefficients
from 40 configurations equilibrated at the operating point (chi, chi2) = (0.5, 0.5), measured by
central differences with common random numbers at four window lengths; and effective-temperature
measurements from 200 starts per condition, with mobility from the drift under a small force on one
species and the conjugate fluctuation from the collective coordinate of that species.

The biological comparison (Section 3.11) uses no new simulation: 184 Chlamydomonas axonemes at 1000
frames per second across eight ATP concentrations, 272,974 frames, from ref. [6].

**Convergence.** Equilibration time grows steeply with system size, and an underequilibrated
large box systematically *understates* the largest cluster — mimicking a finite characteristic
cluster size. Every continuation run raised the largest-cluster fraction: +3.4% (N=4000),
+8.8% (N=9000, χ=1.5), +50.9% (N=16000), +86.8% (N=9000, χ=1), +100.1% (N=22,000). No
scale-selection claim below rests on a run that has not been continued and shown to plateau.
Collapse of the within-run and between-seed scatter is a more reliable convergence indicator
here than a slope fit, because cluster counts fluctuate by up to 30% within a run.

## 3. Results

### 3.1 Two previously proposed knobs do not vary reciprocity

The coupling-strength parameter α̂ **cancels exactly** from the antisymmetric/symmetric ratio
|g_i−g_j|/|g_i+g_j|: both sectors scale linearly in it. Numerically the ratio is identical to
six decimal places across α̂ ∈ {0.003, 0.005, 0.010, 0.015} at every separation tested. Sweeping
α̂ produced edge turnover *falling* from 0.601 to 0.282 — opposite to the prediction it was
meant to test — because α̂ raises overall attraction and condenses the system (n_cl 40 → 9).

The steric size-ratio sweep also fails, because the protocol of ref. [1] pins the EHD radii
while varying steric radii; it therefore varies packing, not reciprocity. Turnover falls 506×
toward the tail-large regime, but by total condensation into a single cluster (n_cl = 1.0,
lcf = 1.00), not by graded suppression of fission.

This is a caution specific to simulation: in experiment a particle's EHD and steric radii are
the same physical size, so the experimental size ratio does tune nonreciprocity.

### 3.2 Nonreciprocity causes arrested coarsening

Paper scale, 3 seeds per level, all other parameters fixed:

| χ | 0 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 | mono ref |
|---|---:|---:|---:|---:|---:|---:|---:|
| n_cl | 11.1 | 2.2 | 17.3 | 56.5 | 95.2 | 137.3 | 11.2 |
| lcf | 0.301 | 0.910 | 0.898 | 0.699 | 0.717 | 0.771 | 0.46 |
| Q | 0.930 | 0.911 | 0.910 | 0.909 | 0.906 | 0.901 | 0.936 |
| σ²_v | 4.19e-10 | 4.28e-10 | 3.65e-09 | 1.37e-08 | 2.76e-08 | 5.57e-08 | 3.07e-10 |

Trend tests across χ (18 runs): n_cl ρ=+0.931 (p=2.1×10⁻⁸); σ²_v ρ=+0.975 (p=7.1×10⁻¹²).
χ = 0 → 1.5: n_cl 11.1 ± 1.3 → 137.3 ± 0.6 (a factor of 12.4); σ²_v × 133. We quote effect sizes
with seed-level scatter rather than p-values here: with three blocked seeds per level, and
within-run fluctuations of up to 30% (Section 2.5), formal significance tests on run-level means
produce implausibly small p-values that reflect the variance-reduction of the blocked design rather
than the evidence. The effects are large enough not to need them. Both contrasts *strengthen* at paper scale relative to pilot scale (n_cl 9.8×,
σ²_v 24×), indicating the pilot runs had not reached steady state.

**The antisymmetric sector therefore alters static structure.** n_cl and lcf are single-snapshot
observables; they move by an order of magnitude under a change confined to the antisymmetric
coupling. The premise that a solenoidal coupling cannot alter stationary structure holds only
when ∇·(ρ_eq **v**) = 0, which generic nonreciprocal couplings do not satisfy.

### 3.3 Newman modularity is degenerate on these networks

Across all 14 conditions run, cluster count spans 40× (1 → 40) while Q spans 14%, and the two
are uncorrelated: Spearman ρ = +0.31, p = 0.28. A 3-cluster state (Q=0.869) and a 33-cluster
state (Q=0.854) are indistinguishable. Louvain's resolution limit merges communities below
~√(2E) ≈ 60 nodes, and the degeneracy survives a 4× increase in system size.

Modularity excess over the degree-preserving null **decreases monotonically with
nonreciprocity** (0.440 at χ=0 → 0.383 at χ=1.5) and is *highest* in the monodisperse
equilibrium reference (0.462). A prediction of sustained modularity excess in the nonreciprocal
case is therefore not merely unsupported; the ordering is inverted.

### 3.4 No circulation in coarse observables, but genuine microscopic irreversibility

The signed area rate was measured in four observable planes at every χ, at both system sizes.
Every value is consistent with zero and none is distinguishable from the χ=0 equilibrium null
(all p ≥ 0.38 at paper scale, where nsnap=200 doubles the statistics).

Entropy production, by contrast, is clearly positive. Paired-bias-subtracted, per particle:

| χ | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 |
|---|---:|---:|---:|---:|---:|
| net EPR | −5.8e-5 | −2.2e-5 | +8.0e-4 | +2.0e-3 | +6.5e-3 |

resolved for χ ≳ 0.75 and buried in bias noise below. The evolving-structure sweep fits χ^3.01,
but structure co-varies with χ there. Holding the configuration fixed and varying χ over a
window too short to relax gives net EPR/χ² flat to **1.27× and 1.61×** on two of three
configurations, confirming the expected quadratic law; the third configuration's bias probe (a
single stochastic realisation) came in 2× high and is instrument noise.

**The dynamics are genuinely irreversible, and that irreversibility does not appear as
circulation in any network observable measured.**

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
**10.7× faster** than thermal ones. This is an amplitude effect, not a scaling one: the
instantaneous drift from the propulsive force scales as N^−0.50, the same exponent as thermal
diffusion, because the pair propulsions add incoherently at every size (coherence C ≡ |Σp|/Σ|p|
satisfies C·√N ≈ const from N=5 to N=80).

Consequently cluster count is **non-monotonic**: n_cl dips to 2.2 at χ=0.25 (below the
reciprocal 11.1) as weak activity completes the arrested phase separation, then rises to 137 as
strong activity drives fission. Fluctuations peak in this crossover region (sd/mean of n_cl:
7.2% at χ=0, 29.7% at χ=0.25, 18.6% at χ=0.5, 7.1% at χ=1.5).

### 3.6 Attributing the reciprocal/nonreciprocal difference

The monodisperse-versus-bidisperse comparison confounds reciprocity with polydispersity. The
χ=0 control separates them (paper scale):

| observable | mono | χ=0 | χ=1 | share attributable to reciprocity |
|---|---:|---:|---:|---:|
| edge turnover | 0.148 | 0.260 | 0.378 | **52%** |
| σ²_v | 3.07e-10 | 4.19e-10 | 2.76e-08 | **99.6%** |

σ²_v is a clean reciprocity probe; edge turnover is about half polydispersity.

### 3.7 Observable reliability

Edge Jaccard turnover fails as a dynamical probe on four counts: 74% of its χ=1 value survives
at χ=0 where the current is provably zero; it correlates with morphology (ρ=+0.795 with n_cl);
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

**Largest cluster ∝ N^1.01**, against 0.00 for a finite characteristic size and 1.00 for phase
separation; lcf is flat at 0.78–0.80 across a fourfold range of system size. The fragment
population is extensive (n_cl ∝ N^0.96) with median cluster size 3 at every box. No interior
peak appears at any size tested.

**The replicate at the source specification.** Composition was corrected to the published
5,000:17,000 (22.7% type-I; our earlier runs used 25%, and type-I particles carry 5.1× the EHD
strength, so the excess biases toward condensation):

Four seeds, each carried forward until its largest-cluster fraction stopped drifting:

| seed | t̂ = 0 – 3.6×10⁵ (**source duration**) | 3.6 – 7.2×10⁵ | 7.2×10⁵ – 1.08×10⁶ | drift in final window |
|---:|---:|---:|---:|---:|
| 1 | 0.378 | **0.811** | — | +1.1% |
| 2 | 0.453 | **0.852** | — | +0.8% |
| 3 | 0.443 | 0.531 | **0.549** | +0.0% |
| 4 | 0.556 | 0.709 | **0.870** | −0.0% |
| mean | **0.458 ± 0.037** | | **0.770 ± 0.075** | |

Seeds 3 and 4 were still drifting after the second window (+10.6% and +12.3% across its halves)
and were therefore extended to a third; all four are converged at the values shown. The paired
increase from the source duration to the converged state is **+0.313 ± 0.073, a 68% rise**
(paired t = +4.3, p = 0.024), and every seed moves in the same direction, by between 24% and 115%.

Two features deserve comment. An earlier two-seed estimate gave 0.415 → 0.831; the four-seed
figures are 0.458 → 0.770 with a wider spread, so the two-seed version somewhat overstated the
magnitude. More interestingly, **the converged state is not unique**: three seeds settle at
0.81–0.87 while seed 3 settles at 0.549 and stays there. That heterogeneity is consistent with the
kinetic-arrest picture of Section 3.5 — a configuration that has separated into two large clusters
which then cannot find each other has no route to a single condensate on any accessible timescale.
The asymptotic largest-cluster fraction is therefore configuration-dependent, which is itself a
statement about arrest rather than about coarsening.

Two results follow. First, **composition does not explain the discrepancy**: at convergence
22.7% type-I gives lcf = 0.831 against 0.786 for 25%, a 5.8% difference rather than the 47%
the underequilibrated comparison suggested. Second, **at the source paper's own simulation
duration this reimplementation reproduces its reported state** — a majority of particles outside
the largest cluster, median cluster size 3 — and the same system continued to twice that
duration phase-separates. Both seeds plateau at lcf ≈ 0.83–0.85 over their final third.

The small-cluster state is therefore a long-lived transient of this model rather than its steady
state. This is both a validation of the reimplementation (it reaches the published state under
the published conditions) and a limitation of the model (that state does not persist).

### 3.9 The selected scale is an event-rate crossover, not a stable cluster size

Two further measurements resolve what the system selects. **Fragments** — all clusters except
the largest — have size statistics independent of system size: mean 5.02, 4.78, 5.34 particles
at N = 4,000, 9,000, 16,000, with an N-independent cutoff (p99 = 33, 35, 26) and a constant
fragment mass fraction of ~16%. A scale is therefore selected.

Its nature follows from resolving merge and split events directly. At the production snapshot
interval this is impossible — a small cluster is typically absorbed within one interval, so
tracking reports condensate size rather than a growth increment. Repeating from equilibrated
configurations at Δt=50 (105,708 cluster-observations, 3 seeds) gives per-cluster rates:

| cluster size | ~2 | ~5 | ~10 | ~19 | ~45 | ~74 |
|---|---:|---:|---:|---:|---:|---:|
| split rate | 0.333 | **0.453** | 0.254 | 0.181 | 0.081 | 0.086 |
| merge rate | 0.319 | 0.312 | **0.334** | 0.283 | 0.269 | 0.145 |

The rates cross at **S\* ≈ 7–8 particles**, with splitting dominant below and merging dominant
above. Splitting therefore dominates below S\* and merging above it. We call this an **event-rate
crossover scale** rather than a critical nucleus: a crossing of per-cluster event *rates* is not
by itself a zero of the size drift ⟨ΔS|S⟩/Δt, and a critical nucleus would properly be
established by that drift changing sign, or by a committor q(S\*) ≈ ½. Both remain to be measured.
The crossing is nonetheless the scale that separates the subcritical population from the
condensate. The fragment population (mean ≈ 5) is the subcritical
vapour, and its N-independence follows because S\* is set by local energetics.

The steady state is thus a condensate coexisting with a subcritical vapour, separated by a
nucleation barrier. No mechanism caps cluster size, which is consistent with §3.8: the largest
cluster grows as N^1.01 and the system phase-separates.

### 3.10 Tier membership: which variational principle applies

The preceding sections establish what nonreciprocity does to this system. We now ask what kind of
variational description it admits. The question is easy to render vacuous — a framework asserting
that "some variational principle applies" forbids nothing — so we treat tier membership as an
experimental matter. We stress that the following is a working taxonomy of *commonly invoked diagnostics*, not a
established classification: the existence of a single effective temperature in particular is a
nonequilibrium construct in the sense of Cugliandolo and Kurchan [11], and is neither necessary
nor sufficient for near-equilibrium response. With that caveat:

| tier | regime | variational object | signature |
|---|---|---|---|
| I | equilibrium | free energy; MaxEnt | EPR = 0; detailed balance; Boltzmann |
| II | linear response | Rayleighian R = Φ̇ + Ψ, Ψ quadratic in rates | fluxes linear in the forces; Onsager reciprocity among conjugate pairs |
| III | far from equilibrium | no general principle | nonlinear flux–force relations; reciprocity fails |

It is worth noting at the outset that least action does not fail for nonreciprocal systems in the
way sometimes supposed. The Onsager–Machlup action S[x] = ∫ (ẋ − μF)²/4D dt is well defined for
*any* drift field, gradient or not. What fails is the equilibrium corollary — that the stationary
distribution is exp(−βU) and that structure is selected by minimising an energy. That corollary,
and not least action itself, is what the solenoidal proposal of Section 1 relied upon.

We tested three tier-II signatures.

**Quadratic dissipation.** Entropy production per particle, measured with the paired-configuration
protocol of Section 2.4 and bias-subtracted against a same-configuration χ = 0 control:

| χ | 1.5 | 2 | 3 | 5 | 8 |
|---|---:|---:|---:|---:|---:|
| EPR | 6.53e−3 | 1.22e−2 | 2.63e−2 | 7.25e−2 | 1.88e−1 |
| EPR/χ² | 2.90e−3 | 3.04e−3 | 2.93e−3 | 2.90e−3 | 2.93e−3 |

EPR = kχ² with k ≈ 2.93 × 10⁻³, flat to **4.8%** across a fourfold range of drive and a fifteenfold
range of entropy production.

We are careful about what this does and does not test. At *frozen* configuration the quadratic
form is close to guaranteed by construction: the antisymmetric force is exactly linear in χ, the
dynamics are overdamped with configuration-independent mobility, and the Stratonovich heat is
bilinear in force and velocity, so the leading term is μ⟨F_a²⟩ ∝ χ² with the cross term vanishing
by the symmetry of the reciprocal sector. The measurement is therefore best read as a **validation
of the entropy-production estimator** — non-trivial given the discretisation problems documented in
Section 5 — rather than as an independent test of linear response. The physically informative
number is the *evolving-structure* exponent, χ^3.01 (Section 3.4), which is not quadratic and which
reflects the co-variation of structure with drive.

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
| 200 | +8.67e−7 | +2.04 | +2.16e−6 | +5.2 |
| 400 | +1.62e−6 | +3.12 | +1.17e−6 | +2.3 |
| 800 | +2.70e−6 | +4.57 | **+1.18e−8** | **+0.0** |

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

We deliberately stop short of calling this Onsager–Casimir reciprocity. That identification would
require a defined microscopic time-reversal operation with the parities of all time-odd control
parameters specified, and it cannot be inferred retrospectively from an observed sign. Three
further caveats apply. The coefficients are measured about a nonequilibrium operating point rather
than about equilibrium; χ and χ₂ are internal coupling constants in the equations of motion, not
thermodynamic affinities, so the bilinear identity above is a definition rather than a derivation;
and because the two channels were constructed with the same pair-co-propulsion structure, the
antisymmetry may follow from that construction symmetry rather than from microreversibility. A
defensible reciprocity test would require either an equilibrium Green–Kubo derivation or an
explicit Fokker–Planck computation of L_ij for this model, together with a demonstration that the
antisymmetry survives when the two channels are made structurally dissimilar. We report the
measurement and its robustness, not a reciprocity relation.

This result is not merely consistent with the preceding one but explains it. Entropy production is
the quadratic form EPR = Σ L_ij X_i X_j, to which only the symmetric part of L contributes; a purely
antisymmetric off-diagonal block dissipates nothing. The cross-coupling between the two
nonreciprocal channels is therefore *reactive* rather than dissipative, which is precisely why the
entropy production reduces to a single-coefficient quadratic. Two measurements made by independent
routes agree on the structure of the response matrix.

The naive parity argument predicts a symmetric relation: both fluxes have the form Σ f · v with f
even under time reversal and v odd, making both fluxes odd and their signatures equal. The
measurement contradicts this robustly. The most economical hypothesis is that the cross-coupling is
*gyroscopic* — a reactive, Magnus- or Coriolis-like coupling in the space of collective
coordinates, which enters a response matrix antisymmetrically without contributing to dissipation.
That is consistent with the vanishing symmetric part, and it would explain why entropy production
reduces to a single-coefficient quadratic. We offer it as a hypothesis requiring analytic
confirmation, not as an established mechanism, and note that the alternative — that the
antisymmetry is an artifact of having built the two channels with identical structure — is equally
live until the dissimilar-channel control is run.

**Effective temperature.** At equilibrium the Einstein relation fixes D/μ = kT for every degree of
freedom; here mobility is 1/s_i and the noise variance σ²/s_i, so D/μ = σ²/2 exactly for both
species. Out of equilibrium the ratio defines an effective temperature, and a single T_eff shared
across degrees of freedom is a tier-II signature. We measured mobility from the drift under a small
force applied to one species and the conjugate fluctuation from the collective coordinate
X = Σ_{i∈species} x_i, with perturbed and unperturbed runs sharing a noise stream. Linear response
was verified separately: T_eff is independent of the perturbation strength down to f = 3 × 10⁻⁶.

At equilibrium the estimator behaves as it must, returning T_eff/kT = 0.98 ± 0.07 for the large
species and 0.94 ± 0.07 for the small, their difference consistent with zero (t = 1.01), and
independent of the measurement window.

Out of equilibrium, however, **the quantity itself ceases to be well defined**, and this is the
substantive finding rather than any particular value. Applying the same plateau test used for the
Onsager coefficients, T_eff at χ = 1.5 approximately doubles between windows of T = 200 and
T = 400 (from 5.48 to 9.78 kT for the large species). The cause is visible in the exponents:

| condition | ⟨ΔX²⟩ scaling | response scaling | resulting T_eff |
|---|---|---|---|
| χ = 0 | t^0.86 | t^0.94 | t^−0.08 — window-independent |
| χ = 1.5 | **t^1.81** | t^0.99 | **t^+0.83 — grows without bound** |

At equilibrium the fluctuation and response exponents match, so their ratio is a constant and
equals kT. At χ = 1.5 the collective coordinate is strongly **superdiffusive** (α ≈ 1.8, approaching
the ballistic value of 2) while the response remains linear in time, so the apparent diffusion
coefficient grows with the observation window and no window-independent effective temperature
exists.

A3 is therefore neither passed nor failed as posed: out of equilibrium, in this system, the
fluctuation–dissipation ratio is not a temperature. The apparent species dependence at χ = 1.5
(t = 3.44 at a single window) is a symptom of the same non-convergence rather than evidence for two
distinct temperatures, and we do not interpret it as such. The physical content is that
nonreciprocity drives the collective coordinate ballistic — directed cluster motion, consistent
with the unjamming mechanism of Section 3.5 — which is precisely the circumstance under which an
effective-temperature description breaks down.

**Interpretation.** The system passes the tier-II tests despite being nonreciprocal, dissipative,
and strongly self-organising. A Rayleighian — minimising Φ̇ + Ψ with Ψ quadratic in the rates — is
therefore the appropriate variational object for its dynamics, and the antisymmetric sector belongs
in the dissipation functional rather than in an action whose extremum is expected to leave the
stationary density invariant.

One boundary should be stated explicitly, because our own data mark it. Entropy production is
quadratic *at fixed structure*, but the structural response to χ is strongly non-linear: cluster
count is non-monotonic, dipping at χ = 0.25 before rising by an order of magnitude (Section 3.5).
The Rayleighian therefore governs the dynamics at fixed structure; *structure selection* is a
non-equilibrium transition that no current variational principle covers. Identifying that gap seems
to us more useful than obscuring it.

### 3.11 A biological comparison: does irreversibility survive coarse-graining?

Section 3.10 has a consequence worth drawing out. If nonreciprocity, dissipation and
self-organisation are jointly insufficient to leave linear response, then they cannot be what
distinguishes living organisation, and this colloid suspension becomes useful as a control: a well
characterised system that is *driven but not alive*. Any proposed signature of biological
organisation must distinguish itself from this system, not merely from equilibrium.

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
a provable equilibrium reference.

We asked whether a biological system differs on that specific axis, using published high-speed
recordings of reactivated *Chlamydomonas* axonemes [6, 22]: 184 usable axonemes imaged at 1000 frames
per second across eight ATP concentrations from 50 to 1000 µM, 272,974 frames in total. Axoneme
shapes were converted to tangent-angle representations, removing rigid translation and rotation,
and reduced by principal component analysis to two dominant shape modes; the signed area rate was
then computed in that plane against phase-randomised surrogates, using the identical estimator
applied to the colloid observables.

We first note that the test originally planned for this dataset — whether entropy production is
quadratic in the ATP drive — is ill-posed. The chemical potential of ATP hydrolysis depends on the
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

Single active colloidal clusters show no circulation either, so the contrast does not stem from the
choice of coarse-graining level. We note the limitation that only three to six clusters per
condition persist long enough to be tracked over the 101-frame window, so the colloid side of this
control rests on few objects, albeit with a null result uniform across all of them.

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
collective order alone nor nonreciprocity alone suffices — both give a null — but a chiral variant,
whose intrinsic turning rate imposes a collective limit cycle, circulates strongly.

Varying the turning rate identifies the origin: the measured area rate is 0.0056, 0.0089, 0.0141
and 0.0352 for ω = 0.005, 0.01, 0.02 and 0.04, a ratio to ω of ≈0.9 throughout. The signal is
therefore largely the *rotation of the mean heading vector at the imposed rate*, not an emergent
property.

We therefore withdraw the interpretation that coarse-grained circulation distinguishes biological
from non-biological organisation. What the measurement detects is whether the system possesses a
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
poor: our strongly nonreciprocal colloid registers none. As a probe of organised, cycle-structured
activity it appears to be informative, because it asks whether directed structure survives
coarse-graining — which is closer to what claims about biological organisation actually assert than
any purely structural measure is.

We are careful about the scope of this claim. Two systems establish that a difference exists on
this axis; they do not establish that the axis separates biological from non-biological systems in
general. Circulation measured in a two-mode projection is moreover a lower bound on irreversibility,
and neither system's full phase space was searched.

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
be large. Nonreciprocity drives arrested coarsening in this model, with a same-particle reciprocal
control that the source experiment cannot itself provide.

### 4.2 The antisymmetric sector restructures rather than merely circulating

The physical result that most changes the interpretation is that the antisymmetric sector alters
static structure, and does so dominantly. A description that assigns the symmetric part to
energy-like structure selection and the antisymmetric part to circulation cannot be maintained:
both sectors determine the structure, and here the antisymmetric one determines most of it.

This had to be so. The premise that a solenoidal coupling preserves the stationary density holds
only when ∇·(ρ_eq **v**) = 0, a condition that generic nonreciprocal couplings do not satisfy. It
also stands in direct tension with the source experiment, whose central finding — that
nonreciprocity arrests coarsening — is itself a report that nonreciprocity changes steady-state
structure.

What the antisymmetric sector does uniquely contribute is a positive quadratic dissipation.
Section 3.10 shows that its cross-coupling to a second nonreciprocal channel is purely reactive, so
that entropy production reduces to a single-coefficient quadratic in the drive. That functional
form is Onsager–Rayleigh structure rather than action-extremum structure, which is why we conclude
that such couplings belong in a dissipation functional.

### 4.3 Network observables are largely blind to nonreciprocity

Several standard network measures fail to detect nonreciprocity here, and they fail in
instructive ways.

Newman modularity is degenerate on these contact graphs. Cluster count varies by a factor of forty
across the pilot-scale conditions while modularity varies by 14%, and the two are uncorrelated; a
three-cluster and a thirty-three-cluster configuration are indistinguishable. The degeneracy
survives a fourfold increase in system size, so it is a property of modularity on two-dimensional
contact networks at this density rather than a finite-size effect. Modularity excess over a
degree-preserving null moreover *decreases* with nonreciprocity and is highest in the reciprocal
reference, inverting the predicted ordering.

Edge turnover fails for four independent reasons: 74% of its value at χ = 1 survives at χ = 0,
where the probability current is provably zero; it correlates strongly with cluster morphology; its
lag-one value is dominated by flicker at the contact threshold rather than by rewiring; and it
becomes non-monotonic in χ at paper scale, a trend reversal that only appeared after
equilibration.

The circulation measure is the most interesting failure, because it turns out to be the right
instrument pointed at the wrong system. It registers nothing in the colloid suspension at either
system size, yet the same estimator applied to cilia gives a median |z| of 3.2. Circulation is
therefore a poor probe of nonreciprocity but a useful probe of organised, cycle-structured
activity.

### 4.4 Two kinds of arrest, distinguished by dissipation

Section 3.5 establishes that the reciprocal reference state is a kinetically arrested gel: pair
wells are 4.8 to 7.0 k_BT deep so bonds do not break thermally, cluster diffusion scales as 1/n so
a thousand-particle cluster moves less than a particle diameter over an entire run, and — decisively
— a condensed configuration evolved at χ = 0 with no activity whatever remains condensed. The
condensate is equilibrium-accessible; the reciprocal dynamics simply cannot reach it.

This carries a caution for the wider literature on this system. Both the reciprocal control and the
monodisperse reference are *arrested*, and they are statistically indistinguishable from one
another. "Arrested coarsening" therefore does not by itself indicate an actively maintained state.
What distinguishes the two regimes is dissipation: exactly zero at χ = 0, positive and quadratic
above it. Where an analogy to homeostatic or biological organisation is intended, the transferable
quantity is an irreversibility measure rather than a structural one, since a frozen aggregate and a
dynamically maintained one may be structurally similar while differing absolutely in dissipation.

### 4.5 Arrest without a steady state

The finite-cluster state of this model is real, reproducible and long-lived — it is what one
observes on the timescale of the source experiment — but it is not asymptotic. Continued to twice
the source duration, the same system phase-separates. A variational account that attributes the
finite-cluster state to a steady-state balance between energy and connection cost is therefore
describing a transient rather than an attractor.

The transient is nonetheless long, and quantifiably so. Fitting the largest cluster as S ~ A t^z
gives z between 0.27 and 0.42 for every nonreciprocal condition, so the crossover time scales as
N^(1/z) with 1/z between 2.4 and 3.7: a tenfold larger system remains in the finite-cluster regime
between 250 and 5000 times longer. In a macroscopic suspension the distinction between long-lived
transient and steady state becomes operationally empty. This reconciles our result with the source
report, which is accurate for its own system size and duration.

A related inversion deserves explicit statement, because the phrase "arrested coarsening" invites
the opposite reading. The reciprocal case has by far the *slowest* coarsening (z = 0.089), and every
nonreciprocal case sits near the diffusion-limited value of 1/3. Nonreciprocity accelerates growth
of the majority phase — it unjams the gel — while sustaining a population of small fragments.
"Arrest" in this model therefore denotes a persistent fragment population, not a slowed majority
phase, and the genuinely arrested state is the reciprocal one.

### 4.6 What the biological comparison does and does not establish

The colloid suspension passes the two tier-II diagnostics that admit a definite answer, while the
third proves ill-posed out of equilibrium, and it does so while being nonreciprocal, dissipative
and strongly self-organising. Those ingredients are therefore jointly insufficient to leave linear
response, and cannot be what distinguishes living organisation.

Cilia differ from the colloid on a measurable axis — their irreversibility survives coarse-graining
— and it is tempting to read that as a signature of biological organisation. A synthetic control
shows that it is not. A chiral Vicsek flock, given a collective limit cycle by construction,
circulates more strongly than cilia do, while a nonreciprocal but non-chiral flock circulates not at
all despite higher polar order. The discriminating property is therefore the presence of a cyclic
collective mode, not biology and not nonreciprocity.

We regard this as the more useful outcome. It converts a claim we could not have defended into a
narrower one that the data support: microscopic irreversibility need not project onto macroscopic
observables, and whether it does is a structural property of the system's collective modes rather
than of its provenance. It also identifies what a genuine biological signature would have to do —
distinguish an emergent limit cycle, such as the coordinated ciliary beat, from an imposed one,
such as a single-particle turning rate — which the area-rate estimator by construction cannot.

### 4.7 Consequences for a network-weighted action principle, and for minimum-action learning

One specific formulation of the decomposition tested above is the Network-Weighted Action Principle
(NWAP), whose directed-graph extension makes the symmetric/antisymmetric assignment explicit [9].
*That framework is the present author's own, and this study was undertaken to test one of its
predictions; we state the outcome directly and note the conflict of interest.* It is mixed, and the parts that fail and the parts that survive are
informative in different ways.

**The structural proposal is validated.** NWAP's directed-graph extension holds that a
nonreciprocal coupling should be decomposed into symmetric and antisymmetric sectors, with the
antisymmetric sector carrying the nonequilibrium content. That decomposition is exactly what this
system required, and constructing it is what made every subsequent measurement possible. Its
central quantitative prediction — that cluster statistics are controlled by the
antisymmetric-to-symmetric ratio — is borne out across two system sizes, with a same-particle
reciprocal control that the source experiment cannot itself provide.

We state this more carefully than a rank correlation would suggest. The Spearman coefficient across
the χ sweep is ρ = +0.931, but the underlying relation is *not monotonic*: cluster count dips to 2.2
at χ = 0.25, a factor of five below its reciprocal value, before rising by two orders of magnitude
(Section 3.5). A rank correlation is a poor summary of such a curve. The defensible statement is
that nonreciprocity strongly controls cluster statistics, that the control is non-monotonic, and
that no proposed functional form — NWAP's included — predicts the dip. This is, to our knowledge, the first
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
size, against a provable equilibrium null. The "circulating partition" half of the prediction has
no empirical support in this system, notwithstanding that the system is genuinely irreversible.

*The modularity signature.* The proposed empirical test — sustained excess Newman modularity in the
nonreciprocal case relative to a reciprocal null — fails with inverted sign. Modularity excess
*decreases* monotonically with nonreciprocity and is highest in the reciprocal reference. The
underlying reason is that Newman modularity is degenerate on these contact graphs, so the test
could not have discriminated in either direction.

**What the antisymmetric sector actually contributes is dissipation.** Section 3.10 establishes
that it produces a positive entropy production quadratic in the coupling, and that its
cross-coupling to a second nonreciprocal channel is purely reactive. That functional form —
a positive quadratic in the drive, with an antisymmetric reactive cross-block — is Onsager–Rayleigh
structure, not action-extremum structure.

**Where this leaves the framework.** We would draw three conclusions, offered constructively.

First, the broadening that NWAP requires is not ad hoc. There is a standard hierarchy of
variational principles for stochastic dynamics, and least action does not in fact fail for
nonreciprocal systems: the Onsager–Machlup action is well defined for any drift field. What fails
is the equilibrium corollary. NWAP's Triple-Action, with its explicit information term traded
against an energy term under structural constraints, has the formal shape of a Maximum Caliber
functional — the path-entropy principle that reduces to MaxEnt statically and to Onsager near
equilibrium. Placing NWAP as a constrained member of that family, with network-structural rather
than thermodynamic constraints, would give it a principled home and two checkable reductions: it
should reduce to a Rayleighian at weak drive and to free-energy minimisation at zero drive.

Second, the discipline that keeps such a broadening falsifiable is that each tier carries its own
experimental signature and membership must be *measured*. We tested three such signatures here.
A framework asserting that "some variational principle applies" forbids nothing; one asserting
tier-II membership predicts quadratic dissipation, Onsager reciprocity, and a single effective
temperature, each of which can fail independently. In this system the first two hold and the third
turns out not to be well defined — an outcome that no amount of reinterpretation could have
produced from the framework alone.

Third, the boundary of applicability should be stated rather than obscured. Entropy production here
is quadratic *at fixed structure*, while the structural response to the antisymmetric coupling is
strongly non-linear. A Rayleighian therefore covers the dynamics at fixed structure; structure
selection is a non-equilibrium transition that no current variational principle covers. Since
structure selection is precisely what NWAP was built to explain, this identifies the gap that new
theory would have to fill, and it is a more useful result for the framework than a claim of
coverage would have been.

**Connection to minimum-action learning.** The same programme has a computational arm, Minimum-Action
Learning (MAL) [9], which selects symbolic force laws from noisy trajectories by minimising a
Triple-Action functional whose discriminating component is an energy-conservation term: among
candidate basis functions, the one that conserves energy under dynamical rollout is selected, and
this criterion is what lifts raw basis-identification rates to complete identification on the Kepler
and Hooke benchmarks.

The present results delimit that method's domain in a specific and, we think, useful way. The
energy-conservation criterion presupposes a conserved energy, which exists precisely when the force
is reciprocal. At χ = 0 the system studied here satisfies that condition exactly; at χ > 0 it does
not, and no amount of basis refinement will recover a conserved quantity that the dynamics do not
possess. MAL applied to trajectories from this model should therefore succeed at χ = 0 and fail
progressively as χ increases — and, more interestingly, the *manner* of its failure is diagnostic:
the conservation residual it computes is, up to normalisation, a measure of the very dissipation we
quantify in Section 3.10, which is quadratic in the antisymmetric coupling at fixed structure.

This suggests a concrete extension rather than merely a limitation. Replacing the
energy-conservation term with a *dissipation-consistency* term — requiring that the selected force
law reproduce the measured entropy production, not that it conserve energy — would carry the method
into nonreciprocal systems while retaining the feature that gives it its discriminating power. The
trajectories and measured dissipation rates generated here constitute a ready benchmark for that
test, and we regard it as the most direct way to connect the two arms of the programme.

**A note on the empirical hook.** The reanalysis originally proposed to connect NWAP to this
literature — computing modularity excess on published trajectory data — should not be attempted;
we ran it and it fails with inverted sign for reasons intrinsic to the measure. The analysis that
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

**Entropy production rests on a subtracted artifact.** The Stratonovich heat estimator subtracts
two nearly cancelling terms of order μ|F|², and its discretisation residual exceeds the signal at
the production timestep. It is usable only at much smaller timestep with a paired
same-configuration control, and it resolves only for χ ≳ 0.75. The bias probe is a single
stochastic realisation and failed on one of three configurations in the pilot measurement.

**Response coefficients required a plateau test.** Both the Onsager and effective-temperature
measurements produced stable, small-error-bar values under protocols that were subsequently shown
to be wrong, and in each case only a scan against the measurement window or the perturbation
strength revealed the problem. We therefore report no response coefficient that has not been shown
to plateau. Readers should treat any such coefficient reported without that check, in this
literature generally, with corresponding caution.

**The analysis graph is undirected.** Contact graphs are constructed with `nx.Graph` throughout,
so the directed-graph content of the motivating proposal is untested by this pipeline. Building the
directed contact graph — with edges signed by the nonreciprocal force imbalance — remains the
natural next step for that specific claim.

**Box geometry in the source-specification replicate.** The source domain is 648 × 360 µm, an
aspect ratio of 1.8, whereas the replicate uses an equal-area square. We tested this directly with
a rectangular-box integrator (`simulate_rect.py`, validated as bitwise identical to the square
kernel in the L_x = L_y limit): at N = 4000, χ = 1.5, matched area, density and duration, the 1.8:1
rectangle gives lcf = 0.809 ± 0.009 against 0.779 ± 0.004 for the square, a difference of 3.8%.
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
indefinitely.

**Sample sizes in the biological comparison.** The cilia analysis uses two shape modes and
phase-randomised surrogates; a stricter surrogate (iterative amplitude-adjusted Fourier transform),
which preserves the amplitude distribution as well as the power spectrum, would make the null more
conservative and would likely reduce the reported z-scores. Beat frequency is a spectral-peak
estimate and is not monotonic at the top of the ATP range.

## 6. Conclusions

1. Isolating the antisymmetric sector of a nonreciprocal coupling requires a parameter that varies
   it alone. Neither the coupling strength nor the steric size ratio does so; the reciprocity
   mixing parameter χ does, and χ = 0 provides a rigorous equilibrium reference rather than merely
   a control.

2. Nonreciprocity causes arrested coarsening in this model: cluster count rises by a factor of 12.4
   at two system sizes with the reciprocal coupling held bit-for-bit fixed.

3. The antisymmetric sector alters *static* structure and is not structure-preserving. The
   solenoidal premise fails, and would have had to, since it contradicts the source experiment's own
   central finding.

4. Newman modularity is degenerate on these networks and its excess over a degree-preserving null
   runs opposite to the predicted direction. Cluster-size statistics, not modularity, carry the
   signature.

5. No circulation is detectable in coarse network observables, yet entropy production is positive
   and quadratic in the antisymmetric coupling. The system is genuinely irreversible, but that
   irreversibility does not survive coarse-graining.

6. The reciprocal reference state is a kinetically arrested gel, not an equilibrium structure. Weak
   nonreciprocity unjams it, producing a non-monotonic dependence of cluster count on χ, and the
   genuinely arrested state is the reciprocal one.

7. No finite characteristic cluster size is selected: the largest cluster grows as N^1.01 across a
   fourfold range of system size. The scale that is selected, at approximately eight particles, is
   a critical nucleus — an unstable fixed point separating a subcritical vapour from an unbounded
   condensate — rather than a stable organisational level. The finite-cluster state is a long-lived
   transient whose lifetime scales as N^(1/z) with 1/z between 2.4 and 3.7.

8. The system occupies the linear-response tier on the two membership tests that admit a definite
   answer: entropy production is quadratic to within 4.8%, and the differential cross-response between two
   nonreciprocal channels is predominantly antisymmetric, with the symmetric part vanishing at
   every measurement window. Whether that constitutes Onsager–Casimir reciprocity is not
   established here. The
   third test is inconclusive for an informative reason — out of equilibrium the collective
   coordinate is superdiffusive (⟨ΔX²⟩ ~ t^1.8), so the fluctuation–dissipation ratio grows with
   the observation window and is not a temperature. A Rayleighian, not an action, is the
   appropriate variational object, and its scope is the dynamics at fixed structure, since
   structure selection itself remains uncovered.

9. Nonreciprocity, dissipation and self-organisation are jointly insufficient to leave linear
   response. This system is consequently a control for systems that are driven but not alive.

10. For the Network-Weighted Action Principle specifically, the structural proposal survives and
    its quantitative prediction is borne out though non-monotonically, while the mechanism attributed to
    the antisymmetric sector — solenoidality, circulation, and a modularity signature — fails on all
    three counts. The sector's measurable contribution is a quadratic dissipation with a reactive
    cross-coupling, which places it in a dissipation functional rather than an action.

11. Cilia show macroscopic circulation in shape space where the colloid, at both system-averaged
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
measured exponent is 1.01. (b) The largest-cluster fraction is flat across a fourfold range of N,
which is the same statement without a fit. `fig3_scaling.png`

**Figure 4 — The event-rate crossover.** (a) Per-cluster split and merge rates against cluster
size, measured at Δt = 50 from equilibrated configurations (105,708 cluster observations). The
rates cross at S* ≈ 7–8. (b) Their difference, showing the sign change. Splitting dominates below
S*, merging above. `fig4_rates.png`

**Figure 5 — Dissipation and response.** (a) Entropy production per particle against χ at fixed
configuration, with the quadratic form. (b) Symmetric and antisymmetric parts of the differential
cross-response between two nonreciprocal channels, against measurement window: the symmetric part
is consistent with zero at every window while the antisymmetric part is stable. (c) Mean-squared
displacement of the collective coordinate against window length: diffusive at χ = 0, strongly
superdiffusive at χ = 1.5, which is why no window-independent effective temperature exists.
`fig5_thermo.png`

**Figure 6 — Circulation at matched level of description.** (a) Median |z| of the signed area rate
for a single tracked colloidal cluster at three values of χ, against a single *Chlamydomonas*
axoneme, using the identical estimator, shape-descriptor reduction to two modes, and
phase-randomised surrogate protocol. (b) The fraction of objects exceeding |z| = 2. Single active
colloidal clusters show no circulation; axonemes show it in 92% of cases. `fig6_circulation.png`

Supplementary figures: `modularity_test.png` (baseline four-panel), `snapshots.png` (final-frame
renders), `coarsening.png` (trajectories by χ), `chi_test.png`, and `sweep.png` (the two
mis-specified sweeps, retained for the record).

## Data and code availability

All code, per-run summaries and figures are in the repository. `simulate.py` (model, reciprocity
parameter, sweeps), `analyze.py` (network observables), `epr_probe.py` (entropy production),
`onsager.py` and `onsager_equil.py` (reciprocity), `fdt.py` (effective temperature),
`arrest_test.py` (configuration-swap test), `crossover.py` (per-cluster propulsive force),
`continue_run.py` (run continuation), `irreversibility.py` and `benchmark_irrev.py` (time-series
estimators and their calibration), `cilia_analysis.py` (biological comparison),
`tests/test_reciprocity.py` (validation of the χ construction).

Supporting documents: `AUDIT.md` (observable validity), `EXPERIMENT.md` (pre-registered protocol
with thresholds committed before running), `BOX_SCALING.md`, `COARSENING.md`, `ONSAGER_RESULT.md`,
`CILIA_RESULT.md`, `VARIATIONAL_FAMILY.md`, `FOLLOWUP.md`, `TODO.md`.

Raw trajectories (approximately 4 GB) are not version-controlled. The cilia dataset is openly
available from Dryad under CC0 [6].

## References

[1] S. Hara, Y. Sumino, et al. *Arrested coarsening in active colloidal suspensions driven by
nonreciprocal electrohydrodynamic interactions.* Phys. Rev. Lett. **137**, 068302 (2026).
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

[9] M. G. Frasch. *Minimum-Action Learning: energy-constrained symbolic model selection for
physical law identification from noisy data.* arXiv:2603.16951 (2026). *(Author's own prior work;
see the disclosure in Section 4.7.)*

[10] E. T. Jaynes. *Macroscopic prediction*, in Complex Systems — Operational Approaches, Springer
(1985), for Maximum Caliber; see also P. D. Dixit et al., J. Chem. Phys. **148**, 010901 (2018).

[11] L. F. Cugliandolo. *The effective temperature.* J. Phys. A **44**, 483001 (2011).

[12] C. P. Broedersz et al. (J. Gladrow, N. Fakhri, F. C. MacKintosh, C. F. Schmidt, D. A. Weitz,
C. Battle, V. F. Geyer, J. Howard). *Broken detailed balance at mesoscopic scales in active
biological systems.* Science **352**, 604 (2016).

[13] F. S. Gnesotto, F. Mura, J. Gladrow, C. P. Broedersz. *Broken detailed balance and
non-equilibrium dynamics in living systems: a review.* Rep. Prog. Phys. **81**, 066601 (2018).

[14] J. Gladrow, C. P. Broedersz, C. F. Schmidt. *Nonequilibrium dynamics of probe filaments in
actin–myosin networks*, and related work on information loss under coarse-graining of
irreversibility. Phys. Rev. E **96**, 022408 (2017).

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
