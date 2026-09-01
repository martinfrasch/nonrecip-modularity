# Which variational principle governs a nonreciprocal active system? Isolating the antisymmetric sector of an electrohydrodynamic colloid model, with a biological comparison

**Martin G. Frasch**

*Computational study, August 2026. Repository: `nonrecip-modularity`.*

---

## Abstract

Nonreciprocal interactions, in which the force particle *i* exerts on *j* is not the negative of
the force *j* exerts on *i*, cannot be derived from a scalar potential. They therefore sit
outside the reach of the variational principles that organise most of equilibrium and
near-equilibrium physics, and it is not obvious a priori which, if any, variational description
survives. The question is sharpened by the Network-Weighted Action Principle (NWAP), a proposed
cross-scale organising principle in which structure is selected by minimising a functional that
trades information against energy under network-structural constraints [9]. NWAP's directed-graph
extension makes a specific, falsifiable prediction about nonreciprocal systems, and this work was
undertaken to test it. We do so in a specific, experimentally grounded system: the
agent-based model of Hara et al. (Phys. Rev. Lett. **137**, 068302, 2026) for size-asymmetric
colloids driven by electrohydrodynamic flows, in which nonreciprocal pair propulsion has been
identified as the minimal ingredient for arrested coarsening.

Our central methodological device is a reciprocity mixing parameter χ that scales the
antisymmetric part of the pair coupling while leaving the symmetric part bit-for-bit unchanged.
The parameter is constructed so that χ = 0 restores Newton's third law exactly and χ = 1 recovers
the published force law; we verify that at χ = 0 the fluctuation–dissipation relation holds, the
stationary state is Boltzmann, and the probability current vanishes identically, so that χ = 0
furnishes a rigorous equilibrium reference rather than merely a convenient control. Two parameters
previously proposed for this role are shown to be unsuitable: the coupling strength α̂ cancels
exactly from the antisymmetric-to-symmetric ratio, and the steric size ratio, as parameterised in
the source model, varies packing while leaving the electrohydrodynamic radii fixed.

Across 129 simulations spanning four system sizes (N = 1,000 to 22,000) we report the following.
First, nonreciprocity is the causal driver of arrested coarsening: cluster count rises from 11.1
to 137.3 (a factor of 12.4, p = 9.4 × 10⁻¹⁰) when the antisymmetric coupling alone is scaled.
Second, and contrary to the assumption that a solenoidal coupling leaves stationary structure
invariant, the antisymmetric sector alters *static* structure by an order of magnitude. Third,
Newman modularity is degenerate on these contact networks: cluster count varies by a factor of 40
while modularity varies by 14%, the two being statistically uncorrelated, and modularity excess
over a degree-preserving null *decreases* with nonreciprocity. Fourth, no circulation is
detectable in any coarse network observable, yet entropy production is positive and scales as χ²,
establishing genuine irreversibility whose signature does not survive coarse-graining. Fifth, the
reciprocal reference state proves to be a kinetically arrested gel rather than an equilibrium
structure, and weak nonreciprocity unjams it, producing a non-monotonic dependence of cluster
count on χ. Sixth, no finite characteristic cluster size is selected: the largest cluster grows in
proportion to system size, and the scale that the system does select — approximately eight
particles — is a critical nucleus, an unstable fixed point, rather than a stable attractor.

We then place the system within the standard hierarchy of variational principles by measuring
tier-membership signatures directly. Entropy production is quadratic in the antisymmetric coupling
to within 4.8% across a fourfold range of drive; the Onsager cross-coefficients between two
independent nonreciprocal channels satisfy reciprocity in its antisymmetric (Casimir) form, with
the symmetric part of the cross-coupling vanishing at every measurement window; and the
fluctuation–dissipation ratio returns exactly kT at equilibrium. The third test is inconclusive out
of equilibrium for a specific and informative reason: the collective coordinate becomes
superdiffusive, so the fluctuation–dissipation ratio grows with the observation window and is not a
temperature at all. The system is therefore in the linear-response tier on the two tests that admit
a definite answer, despite being nonreciprocal, dissipative, and strongly self-organising.

That last finding motivates a comparison with a biological system. Using published high-speed
recordings of reactivated *Chlamydomonas* axonemes across eight ATP concentrations, we show that
cilia differ from the colloid suspension on precisely the axis where the colloid is silent: their
irreversibility *survives coarse-graining*, producing macroscopic circulation in shape space
(median |z| = 3.2 over 184 axonemes, with 92% exceeding |z| = 2) where the colloid shows none. We
argue that the discriminating property of biological organisation, at least in this comparison, is
not directed influence per se — the colloid has ample directed influence — but directed influence
that persists under coarse-graining.

**Keywords:** nonreciprocal interactions, active matter, arrested coarsening, entropy production,
Onsager reciprocity, effective temperature, kinetic arrest, variational principles

---
## 1. Introduction

The principle of least action organises much of physics, but it presupposes that forces derive
from a potential. Nonreciprocal interactions violate that premise by construction. When the force
that *i* exerts on *j* is not the negative of the force *j* exerts on *i*, the force field is not
a gradient, no scalar potential exists, and the equilibrium apparatus that follows from one — a
Boltzmann stationary distribution, structure selected by energy minimisation, detailed balance —
does not apply. Nonreciprocal systems are consequently a natural place to ask which variational
description, if any, survives.

The question has become concrete rather than merely formal. Hara, Sumino and colleagues recently
reported a controlled experimental realisation: polystyrene colloids of two radii, confined
between indium tin oxide electrodes under an AC field, develop electrohydrodynamic flows whose
strength scales steeply with particle radius [1]. Size-asymmetric pairs therefore experience
imbalanced attraction and spontaneously self-propel. The resulting suspension exhibits *arrested
coarsening*: clusters continuously fragment and reorganise rather than growing without bound as
they do in the monodisperse case. Accompanying agent-based simulations identify nonreciprocal pair
propulsion as the minimal ingredient for this behaviour.

A specific proposal motivated the present work: the **Network-Weighted Action Principle** (NWAP),
developed as a candidate cross-scale organising principle in which structure is selected by
extremising a functional that trades information against energy subject to network-structural
constraints [9]. In its implemented form the functional is a "Triple-Action" combining an
information term I_max, an energy term E_min, and a symmetry constraint; schematically
S_NW = ∫ (E − I + A·C) dt, with C a connection-cost term. The framework's stated ambition is to
apply across scales, and it has been argued to generate modular architectures of the kind seen in
evolved systems.

NWAP as originally formulated places its weights on an *undirected* network, which on variation
yields reciprocal forces and therefore cannot describe the present system at all. The proposed
remedy — and the object of this study — is a directed-graph extension in which the adjacency G_ij
is decomposed into symmetric and antisymmetric parts. The symmetric part is held to generate the
gradient, energy-like component of the dynamics; the antisymmetric part is held to generate a
*solenoidal*, circulating component, identified with exactly the nonreciprocal propulsion that
Hara et al. isolate as their minimal ingredient.

Three consequences follow from that assignment, and all three are testable. First, because a
solenoidal field is divergence-free, the antisymmetric sector should not alter any static
structural observable. Second, its entire signature should therefore appear in probability
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
temperature (Section 3.10). The system passes all three despite being nonreciprocal, dissipative
and strongly self-organising, which places it in the linear-response tier and identifies a
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

### 2.3 χ=0 is an exact equilibrium reference

At χ=0 all forces are central, pairwise and equal-and-opposite, hence conservative; mobility
μ_i = 1/s_i with noise variance σ²/s_i gives D_i/μ_i = σ²/2, uniform across particles.
Fluctuation–dissipation therefore holds, detailed balance holds, and **the probability current
is exactly zero**. This is a stronger reference than the monodisperse comparison used
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
χ=0 → 1.5: n_cl 11.1 → 137.3 (12.4×, Welch t=−302.6, p=9.4×10⁻¹⁰); σ²_v ×133
(p=2.9×10⁻⁵). Both contrasts *strengthen* at paper scale relative to pilot scale (n_cl 9.8×,
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

| window | lcf | seed spread | within-run sd | n_cl |
|---|---:|---:|---:|---:|
| t̂ = 0 – 3.6×10⁵ (**the source paper's duration**) | **0.415** | 0.038 | 0.125 | 431 |
| t̂ = 3.6×10⁵ – 7.2×10⁵ | **0.831** | 0.021 | 0.019 | 458 |

Two results follow. First, **composition does not explain the discrepancy**: at convergence
22.7% type-I gives lcf = 0.831 against 0.786 for 25%, a 5.8% difference rather than the 47%
the underequilibrated comparison suggested. Second, **at the source paper's own simulation
duration this reimplementation reproduces its reported state** — a majority of particles outside
the largest cluster, median cluster size 3 — and the same system continued to twice that
duration phase-separates. Both seeds plateau at lcf ≈ 0.83–0.85 over their final third.

The small-cluster state is therefore a long-lived transient of this model rather than its steady
state. This is both a validation of the reimplementation (it reaches the published state under
the published conditions) and a limitation of the model (that state does not persist).

### 3.9 The selected scale is a critical nucleus, not a stable cluster size

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
above. This is an *unstable* fixed point — a critical nucleus — so clusters below S\* dissolve
and clusters above it grow without bound. The fragment population (mean ≈ 5) is the subcritical
vapour, and its N-independence follows because S\* is set by local energetics.

The steady state is thus a condensate coexisting with a subcritical vapour, separated by a
nucleation barrier. No mechanism caps cluster size, which is consistent with §3.8: the largest
cluster grows as N^1.01 and the system phase-separates.

### 3.10 Tier membership: which variational principle applies

The preceding sections establish what nonreciprocity does to this system. We now ask what kind of
variational description it admits. The question is easy to render vacuous — a framework asserting
that "some variational principle applies" forbids nothing — so we treat tier membership as an
experimental matter. The standard hierarchy for stochastic dynamics carries distinct, falsifiable
signatures at each level:

| tier | regime | variational object | signature |
|---|---|---|---|
| I | equilibrium | free energy; MaxEnt | EPR = 0; detailed balance; Boltzmann |
| II | linear response | Rayleighian R = Φ̇ + Ψ, Ψ quadratic in rates | EPR ∝ (drive)²; Onsager reciprocity; single effective temperature |
| III | far from equilibrium | no general principle | EPR non-quadratic; reciprocity fails; observable-dependent temperature |

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
range of entropy production. The relation could have failed at four independent points and did not.

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
part is large and stable. That is, L₁₂ = −L₂₁: reciprocity holds in the Onsager–Casimir
antisymmetric form, the form appropriate when the conjugate variables carry opposite time-reversal
signature.

This result is not merely consistent with the preceding one but explains it. Entropy production is
the quadratic form EPR = Σ L_ij X_i X_j, to which only the symmetric part of L contributes; a purely
antisymmetric off-diagonal block dissipates nothing. The cross-coupling between the two
nonreciprocal channels is therefore *reactive* rather than dissipative, which is precisely why the
entropy production reduces to a single-coefficient quadratic. Two measurements made by independent
routes agree on the structure of the response matrix.

We note one unresolved point. The naive parity argument predicts a symmetric relation: both fluxes
have the form Σ f · v with f a function of positions (even under time reversal) and v odd, making
both fluxes odd and their signatures equal. The measurement robustly contradicts this across three
window lengths. Either the effective time-reversal signature of the two channels differs for a
reason the naive argument misses, or the coupling has a geometric origin outside the standard
parity assignment. We record this as an open theoretical question rather than resolving it here.

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

We asked whether a biological system differs on that specific axis, using published high-speed
recordings of reactivated *Chlamydomonas* axonemes [6]: 184 usable axonemes imaged at 1000 frames
per second across eight ATP concentrations from 50 to 1000 µM, 272,974 frames in total. Axoneme
shapes were converted to tangent-angle representations, removing rigid translation and rotation,
and reduced by principal component analysis to two dominant shape modes; the signed area rate was
then computed in that plane against phase-randomised surrogates, using the identical estimator
applied to the colloid observables.

We first note that the test originally planned for this dataset — whether entropy production is
quadratic in the ATP drive — is ill-posed. The chemical potential of ATP hydrolysis is
approximately 20 k_BT at every concentration in the series, so cilia are far from equilibrium by
construction rather than by measurement, and the linear-response question answers itself.

| [ATP] µM | n | beat frequency (Hz) | median \|z\| | fraction \|z\|>2 | \|area per cycle\| |
|---:|---:|---:|---:|---:|---:|
| 50 | 10 | 14.5 | 3.3 | 1.00 | 6.07 |
| 100 | 11 | 27.7 | 3.0 | 0.91 | 6.23 |
| 240 | 19 | 46.3 | 2.7 | 0.89 | 6.15 |
| 500 | 13 | 91.7 | 4.5 | 0.92 | 5.95 |
| 1000 | 43 | 65.0 | 3.3 | 0.91 | 6.09 |

Pooled over all 184 axonemes, the median |z| is 3.2, with 92% exceeding |z| = 2 and 55% exceeding
|z| = 3. The colloid control, using the same estimator on coarse observables, gives |z| < 1.5 with
nothing above the null. Both systems are irreversible microscopically; only the biological one has
irreversibility that survives coarse-graining into macroscopic circulation.

A second and unanticipated result emerges from the same analysis. The area enclosed per beat cycle
is 5.71 to 6.24 in standardised shape coordinates at every ATP concentration — flat across a
twentyfold range of concentration and a sixfold range of beat frequency — and that value is
approximately 2π, the geometric maximum for a circular orbit in these coordinates. The beat
therefore traces a nearly perfect limit cycle whose *shape* is saturated and ATP-independent, while
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
across the conditions we ran while modularity varies by 14%, and the two are uncorrelated; a
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

The colloid suspension passes all three tier-II tests while being nonreciprocal, dissipative and
strongly self-organising. Those ingredients are therefore jointly insufficient to leave linear
response, and cannot be what distinguishes living organisation.

Cilia differ from the colloid on a specific and measurable axis: their irreversibility survives
coarse-graining. This suggests that the discriminating property is not directed influence as such —
the colloid has ample directed influence, with provably positive entropy production — but directed
influence that persists to macroscopic scales. If that generalises, the appropriate measurement for
claims about biological organisation is whether directed structure survives coarse-graining, rather
than any static or purely microscopic quantity.

We do not claim more than the comparison supports. Two systems demonstrate that a difference exists
on this axis. Establishing that the axis separates biological from non-biological systems in general
would require many more systems on both sides, chosen so that the biological ones are not simply
more strongly driven.

### 4.7 Consequences for the Network-Weighted Action Principle

Because this study was undertaken to test a specific NWAP prediction, we state the outcome for
that framework directly. It is mixed, and the parts that fail and the parts that survive are
informative in different ways.

**The structural proposal is validated.** NWAP's directed-graph extension holds that a
nonreciprocal coupling should be decomposed into symmetric and antisymmetric sectors, with the
antisymmetric sector carrying the nonequilibrium content. That decomposition is exactly what this
system required, and constructing it is what made every subsequent measurement possible. Its
central quantitative prediction — that cluster statistics are set by the antisymmetric-to-symmetric
ratio — holds at ρ = +0.931, p = 2 × 10⁻⁸, across two system sizes, with a same-particle reciprocal
control that the source experiment cannot itself provide. This is, to our knowledge, the first
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
aspect ratio of 1.8; our kernel assumes a square box, so the replicate uses an equal-area square.
A condensate growing along the shorter dimension of the true domain could behave differently. This
is a deviation, not a correction.

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
   (p = 9.4 × 10⁻¹⁰) at two system sizes with the reciprocal coupling held bit-for-bit fixed.

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
   answer: entropy production is quadratic to within 4.8%, and Onsager reciprocity holds in its
   antisymmetric (Casimir) form with the symmetric cross-coupling vanishing at every window. The
   third test is inconclusive for an informative reason — out of equilibrium the collective
   coordinate is superdiffusive (⟨ΔX²⟩ ~ t^1.8), so the fluctuation–dissipation ratio grows with
   the observation window and is not a temperature. A Rayleighian, not an action, is the
   appropriate variational object, and its scope is the dynamics at fixed structure, since
   structure selection itself remains uncovered.

9. Nonreciprocity, dissipation and self-organisation are jointly insufficient to leave linear
   response. This system is consequently a control for systems that are driven but not alive.

10. For the Network-Weighted Action Principle specifically, the structural proposal survives and
    its quantitative prediction holds (ρ = +0.931, p = 2 × 10⁻⁸), while the mechanism attributed to
    the antisymmetric sector — solenoidality, circulation, and a modularity signature — fails on all
    three counts. The sector's measurable contribution is a quadratic dissipation with a reactive
    cross-coupling, which places it in a dissipation functional rather than an action.

11. Cilia differ from the colloid on a specific measurable axis: their irreversibility survives
    coarse-graining, producing macroscopic circulation in shape space where the colloid produces
    none. The candidate signature of biological organisation is therefore directed influence that
    persists under coarse-graining, rather than directed influence as such.

## Figures

| file | content |
|---|---|
| `chi_test.png` | turnover and modularity against χ; relative change showing modularity flat while currents rise |
| `coarsening.png` | largest-cluster fraction and cluster count trajectories by χ; the non-monotonic steady state |
| `modularity_test.png` | baseline four-panel comparison, seed-averaged |
| `snapshots.png` | final-frame renders, monodisperse against bidisperse |
| `sweep.png` | the two mis-specified sweeps, retained for the record |

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

[1] K. Hara, Y. Sumino, et al. *Arrested coarsening in active colloidal suspensions driven by
nonreciprocal electrohydrodynamic interactions.* Phys. Rev. Lett. **137**, 068302 (2026).
DOI 10.1103/96ky-d1p9; arXiv:2509.23164.

[2] M. Fruchart, R. Hanai, P. B. Littlewood, V. Vitelli. *Non-reciprocal phase transitions.*
Nature **592**, 363 (2021).

[3] M. E. J. Newman. *Modularity and community structure in networks.* PNAS **103**, 8577 (2006).

[4] S. Fortunato, M. Barthélemy. *Resolution limit in community detection.* PNAS **104**, 36 (2007).

[5] K. Sekimoto. *Stochastic Energetics.* Lecture Notes in Physics 799, Springer (2010).

[6] V. F. Geyer, J. Howard, P. Sartori. *Ciliary beating patterns map onto a low-dimensional
behavioural space.* Dryad dataset, DOI 10.5061/dryad.0gb5mkm2j (2022). Associated publication:
Nature Physics **18**, 1 (2022), DOI 10.1038/s41567-021-01446-2.

[7] L. Onsager. *Reciprocal relations in irreversible processes.* Phys. Rev. **37**, 405 (1931).

[8] H. B. G. Casimir. *On Onsager's principle of microscopic reversibility.* Rev. Mod. Phys. **17**,
343 (1945).

[9] M. G. Frasch. *Minimum-Action Learning: energy-constrained symbolic model selection for
physical law identification from noisy data.* Preprint (2026). Framework materials and the
Network-Weighted Action Principle: minAction.net.

[10] E. T. Jaynes. *The minimum entropy production principle.* Annu. Rev. Phys. Chem. **31**, 579
(1980); and *Macroscopic prediction*, in Complex Systems — Operational Approaches, Springer (1985),
for Maximum Caliber.
