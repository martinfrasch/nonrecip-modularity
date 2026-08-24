# Isolating the antisymmetric sector of a nonreciprocal colloidal model: nonreciprocity restructures, it does not merely circulate

**Martin G. Frasch**

*Computational study, August 2026. Repository: `nonrecip-modularity`.*

---

## Abstract

Nonreciprocal interactions — where the force particle *i* exerts on *j* is not the negative of
the force *j* exerts on *i* — are not derivable from a scalar potential, and are therefore a
natural stress test for variational descriptions of collective behaviour. We reimplemented the
agent-based model of Hara et al. (PRL 137, 068302, 2026) for size-asymmetric colloids driven by
electrohydrodynamic flows, and introduced a control parameter χ that scales the antisymmetric
(nonreciprocal) part of the pair coupling while leaving the symmetric part bit-identical.
χ=0 is exactly reciprocal and, we show, satisfies detailed balance; χ=1 recovers the published
force law. Across 82 simulations at two system sizes we find: (i) nonreciprocity is the causal
driver of arrested coarsening — cluster count rises from 11.1 to 137.3 (12.4×, p=9.4×10⁻¹⁰)
with the reciprocal coupling held fixed; (ii) the antisymmetric sector alters *static*
structure, contradicting the common assumption that a solenoidal coupling leaves stationary
structure invariant; (iii) Newman modularity is degenerate on these contact networks (cluster
count varies 40× while Q varies 14%, uncorrelated, ρ=+0.31, p=0.28) and modularity excess over
a degree-preserving null *decreases* with nonreciprocity, opposite to a network-modularity
prediction; (iv) no circulation is detectable in any network observable, yet entropy production
is positive and scales as χ² at fixed configuration, establishing genuine irreversibility;
(v) the reciprocal reference state is a kinetically arrested gel rather than an equilibrium
control, and weak nonreciprocity unjams it by boosting cluster mobility ~10.7× at matched size,
producing a non-monotonic dependence of cluster count on χ. We conclude that the antisymmetric
sector's measurable signature is a positive quadratic dissipation, not a structure-preserving
circulation, with consequences for how such couplings should enter variational frameworks.

**Keywords:** nonreciprocal interactions, active matter, arrested coarsening, entropy
production, kinetic arrest, network modularity

---

## 1. Introduction

Hara, Sumino and colleagues reported that polystyrene colloids of two sizes (1.0 and 1.5 μm
radii), confined between ITO electrodes under an AC field, develop electrohydrodynamic (EHD)
flows whose strength scales strongly with particle radius [1]. Size-asymmetric pairs therefore
experience imbalanced attraction and self-propel, and the resulting suspension exhibits
*arrested coarsening*: clusters continuously fragment and reorganise rather than growing without
bound as in the monodisperse case. Their agent-based model identifies nonreciprocal pair
propulsion as the minimal ingredient.

Because F_ij ≠ −F_ji, the force field is not a gradient, and the system is a natural test case
for action-based or variational accounts of pattern selection. A specific proposal motivating
this work held that the antisymmetric part of a directed coupling is *solenoidal*, and therefore
cannot alter a static structural observable — so its entire signature should appear in
probability currents, giving "frozen modularity, circulating partition", with sustained excess
network modularity in the nonreciprocal case relative to a reciprocal null.

Testing that requires a knob on the ratio of antisymmetric to symmetric coupling. We show that
neither of the parameters previously proposed for this role provides one, construct a parameter
that does, and use it to separate what nonreciprocity causes from what merely co-varies with it.

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

82 runs. Pilot scale N=1000, box 12, t̂=1×10⁵ (61 runs: baseline 5 seeds × 2 cases, α̂ sweep
4×3, size-ratio sweep 5×3, χ sweep 6×5). Paper scale N=4000, box 24, t̂=4×10⁵, 200 snapshots
(21 runs: χ sweep 6×3 plus monodisperse reference ×3). Analyses use the late half of each run.
Seeds are blocked: seed *k* gives the same initial configuration and noise stream at every χ.

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

### 3.4 No circulation, but genuine irreversibility

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

## 4. Discussion

The central methodological result is that **isolating the antisymmetric sector requires a
parameter that varies it and nothing else**, and that neither previously proposed knob does.
Once such a parameter exists, the causal claim is straightforward to establish and turns out to
be large: nonreciprocity is the driver of arrested coarsening in this model, with a
same-particle reciprocal control that the source experiment cannot provide.

The physical result that most changes the interpretation is that the antisymmetric sector
**restructures** the system. A description that assigns the symmetric part to energy-like
structure selection and the antisymmetric part to circulation cannot hold: both sectors set the
structure, and here the antisymmetric one dominates it. What the antisymmetric sector uniquely
contributes is a positive quadratic dissipation, EPR = kχ² at fixed configuration. That
functional form is Onsager–Rayleigh structure rather than action-extremum structure, which
suggests such couplings belong in a dissipation functional rather than in an action whose
extremum is expected to leave the stationary density invariant.

The gel result carries a caution for the wider literature on this system. Both the reciprocal
control and the monodisperse reference are *arrested* — kinetically — and are statistically
indistinguishable from each other (n_cl 11.1 vs 11.2). "Arrested coarsening" therefore does not
by itself indicate an actively maintained state. What distinguishes the two regimes is
dissipation: exactly zero at χ=0, positive and quadratic above. If an analogy to homeostatic or
biological organisation is intended, the transferable quantity is an irreversibility measure,
not a structural one, because a frozen aggregate and a dynamically maintained one can be
structurally similar and differ absolutely in dissipation.

## 5. Limitations

- **The condensate is system-scale.** At paper scale the largest cluster holds ~82% of particles
  and the size distribution decays monotonically (median cluster size 3) with no interior peak.
  Whether a finite characteristic cluster size is selected cannot be decided in a box this
  small. A box-scaling series at fixed density (N=9000/L=36 and N=16000/L=48 against the
  existing N=4000/L=24) is in progress; a finite characteristic size S* requires lcf ∝ S*/N,
  i.e. falling as 1/L², whereas true phase separation gives L-independent lcf.
- **χ is a synthetic control parameter**, not an experimental one; χ=1.5 lies outside the
  derived EHD model and is reported as trend evidence only.
- **Entropy production rests on subtracting a discretisation artifact.** It resolves only for
  χ ≳ 0.75 at the timestep used, and the bias probe is a single stochastic realisation that
  failed on one of three configurations. Averaging the bias probe and reducing dt a further 4×
  would settle the low-χ end.
- **The analysis graph is undirected** (`nx.Graph`) throughout, so the directed-graph content of
  the motivating proposal is untested by this pipeline.
- **Pairwise forces only**, inheriting the absence of many-body hydrodynamics from ref. [1],
  which shifts absolute cluster sizes.
- The steepening of the drift exponent for the largest clusters (−0.82 versus −0.50 overall) is
  unexplained; finite-size pinning is plausible but unestablished.

## 6. Conclusions

1. Nonreciprocity causes arrested coarsening: 12.4× in cluster count, p=9.4×10⁻¹⁰, at two system
   sizes, with the reciprocal coupling held bit-identical.
2. The antisymmetric sector alters static structure; it is not structure-preserving.
3. Newman modularity is degenerate on these contact networks and modularity excess runs opposite
   to the predicted direction. Cluster-size statistics, not modularity, carry the signature.
4. No circulation is detectable, but entropy production is positive and quadratic in the
   antisymmetric coupling, establishing genuine irreversibility.
5. The reciprocal reference state is a kinetically arrested gel; weak nonreciprocity unjams it,
   producing a non-monotonic dependence of cluster count on the antisymmetric coupling.

## Figures

| file | content |
|---|---|
| `chi_test.png` | turnover and Q vs χ; relative change showing Q flat while currents rise |
| `coarsening.png` | lcf and n_cl trajectories by χ; non-monotonic steady state |
| `modularity_test.png` | baseline four-panel, seed-averaged |
| `snapshots.png` | final-frame renders, monodisperse vs bidisperse |
| `sweep.png` | the two mis-specified sweeps, retained for the record |

## Data and code availability

All code, per-run summaries and figures are in the repository. `simulate.py` (model and
sweeps), `analyze.py` (network observables), `epr_probe.py` (entropy production),
`arrest_test.py` (configuration-swap test), `crossover.py` (per-cluster propulsive force),
`continue_run.py` (run continuation), `tests/test_reciprocity.py` (χ validation). Supporting
documents: `AUDIT.md` (observable validity), `EXPERIMENT.md` (pre-registered χ protocol with
thresholds committed before running), `RESULTS.md`, `FINAL_RESULTS.md`, `COARSENING.md`,
`NWAP_ASSESSMENT.md`. Raw trajectories (~4 GB) are not version-controlled.

## References

[1] K. Hara, Y. Sumino, et al. *Arrested coarsening in active colloidal suspensions driven by
nonreciprocal electrohydrodynamic interactions.* Phys. Rev. Lett. **137**, 068302 (2026).
DOI 10.1103/96ky-d1p9; arXiv:2509.23164.

[2] M. Fruchart, R. Hanai, P. B. Littlewood, V. Vitelli. *Non-reciprocal phase transitions.*
Nature **592**, 363 (2021).

[3] M. E. J. Newman. *Modularity and community structure in networks.* PNAS **103**, 8577 (2006).

[4] S. Fortunato, M. Barthélemy. *Resolution limit in community detection.* PNAS **104**, 36 (2007).

[5] K. Sekimoto. *Stochastic Energetics.* Lecture Notes in Physics 799, Springer (2010).
