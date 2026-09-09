# Plan: cracking the questions the manuscript leaves open

*Drafted 2026-09-08 against commit 074e011. Companion scratch code: `static_epr.py` (validated below).*

The manuscript closes with twelve conclusions and a limitations section that, read together,
leave nine questions genuinely open. This plan takes each one, asks first whether reasoning alone
settles it, then whether existing data settle it, and only then what new data would look like.
The items are ordered by expected yield per unit of compute. Total new simulation for tracks A–G is
under ~15 CPU-hours on the M2 Max; only track H needs the M5 Studio or a different model.

## 0. The open questions, as the manuscript states them

| # | question | where it is left | status in paper |
|---|---|---|---|
| Q1 | Is the two-term dissipation aχ² + bχ with b < 0 real, or a system-size confound? | §3.10, audit §Q, abandoned N = 4000 run | conservative wording, fit range confounded with N |
| Q2 | What does a *negative* linear term mean physically? | §3.10 | "sign matches the construction's cross-term" — no mechanism |
| Q3 | Is there any coordinate on which an effective temperature exists? | §3.10 A3, Concl. 8 | "not a well-defined quantity on this coordinate" |
| Q4 | Is S* ≈ 7–8 a critical nucleus? | §3.9, Concl. 7 | drift ⟨ΔS\|S⟩ and committor "remain to be measured"; no error bar on S* |
| Q5 | What selects structure, i.e. what predicts the dip in n_cl at χ = 0.25? | §3.5, §4.7 | "no proposed functional form predicts the dip" |
| Q6 | Is the fragmented state transient or arrested in the thermodynamic limit? | §3.8, Concl. 7 | lifetime scaling "an extrapolation, not a measured lifetime" |
| Q7 | Does irreversibility project onto *any* coarse observable? | §3.4, Concl. 5, 12 | four planes null; "a statement about those projections" |
| Q8 | The directed contact graph, never built | §5 Limitations | "the natural next step" |
| Q9 | A defensible Onsager / linear-response test | §3.10 | "would require a Green–Kubo derivation … we do not have one" |

Three smaller limitations (σ²_v definition, iAAFT surrogates for cilia, rectangular box at
N = 22 000) are picked up in track G.

## Track A — Dissipation from configurations alone (Q1, Q2; existing data)

**Reasoning result.** The Stratonovich heat is Q̇ = Σ_i F_i∘ẋ_i with ẋ_i = F_i/s_i + ξ_i. Writing
F = F_s + χF_a and using the Itô–Stratonovich conversion ⟨g∘ξ⟩ = D∇·g,

    ⟨Q̇⟩ = ⟨F_s∘ẋ⟩ + χ ⟨F_a∘ẋ⟩,   ⟨F_s∘ẋ⟩ = −d⟨U⟩/dt,
    ⟨F_a∘ẋ⟩ = Σ_i (1/s_i) ⟨F_a,i · F_i⟩ + Σ_i D_i ⟨∇_i · F_a,i⟩,   D_i = σ²/(2 s_i).

Every term on the right is an analytic function of the configuration: F_a is the known
pair-propulsion field, F is the full force, and ∇·F_a is (2h + r h′)(D_i − D_j) per pair plus a
cutoff delta at r = 1. In a stationary (or slowly coarsening, measurably so) state the symmetric
sector drops out entirely. **The excess dissipation is therefore a static ensemble average over
stored snapshots, with no trajectory, no small timestep, no paired bias probe, and no subtraction
of two O(μ|F|²) terms.** This is the mean-drift form of the entropy production; what makes it
usable here is that the non-conservative part of the drift is known exactly.

Two corollaries settle Q2 by reasoning:

1. Split ⟨F_a∘ẋ⟩ = b̃(χ) + χ ã(χ) with ã = Σ⟨|F_a|²⟩/s_i ≥ 0 and b̃ = Σ⟨F_a·F_s⟩/s_i + Σ D⟨∇·F_a⟩.
   In a Boltzmann ensemble b̃ ≡ 0 by integration by parts (μ⟨F_a·(−∇U)⟩_eq = −D⟨∇·F_a⟩_eq). So the
   "linear term" is nonzero only because the χ-ensemble is not Boltzmann: **b̃ measures the
   structural response to the drive**, and its sign is fixed by Le Chatelier-type reasoning to be
   negative (the structure rearranges to absorb the push). It is a quadratic-order effect that the
   fit range makes look linear.
2. EPR ≥ 0 for every χ with EPR(0) = 0 forces EPR′(0) ≥ 0. A genuinely linear term with b < 0 at
   the origin is therefore impossible; the fitted aχ² + bχ goes negative below χ = 0.45 and is the
   wrong parameterisation at low drive. The right statement is EPR = χ² A(χ) with A(χ) = ã + b̃/χ
   rising with χ because unjamming (§3.5) changes the structure.

**Validation already done (this session).** On the χ = 2, 3, 5, 8 runs at N = 4000, late half,
20 snapshots, seed 1, at the production timestep:

| χ | paired probe (paper) | configurational | ratio |
|---|---:|---:|---:|
| 2 | 1.22e−2 | 1.37e−2 ± 0.19e−2 | 1.12 |
| 3 | 2.63e−2 | 3.05e−2 ± 0.56e−2 | 1.16 |
| 5 | 7.25e−2 | 5.63e−2 ± 1.0e−2 | 0.78 |
| 8 | 1.88e−1 | 1.81e−1 ± 0.20e−1 | 0.96 |

Agreement within the (large, single-seed, 20-snapshot) error at every point. Two things the
decomposition shows immediately. First, the |F_a|² term and the F_a·F_s cross-term cancel to ~90%
(e.g. +0.104 and −0.094 at χ = 2), so **about nine tenths of the nonreciprocal push is absorbed by
contact forces and only the residue dissipates** — that is the physical content of the negative
cross-term. Second, the symmetric-sector identity Q̇_s = −dU/dt, which must hold in stationarity,
fails on dt = 0.05 snapshots (Q̇_s ≈ +4.7 T per particle per unit time at χ = 0, against a measured
−dU/dt of 1.5 × 10⁻⁶): the production timestep distorts contact-pair statistics
(dt·k·μ_max = 0.5). Re-equilibrating each snapshot for t = 50 at smaller dt brings the residual
down as 4.7 → 1.1 → 0.45 for dt = 0.05 → 0.0125 → 0.0025, i.e. O(dt) as expected, so the
estimator is sound and needs one calibration per condition, not a paired probe. At dt = 0.005 the
χ = 2 and χ = 8 excess values move by less than their own error (0.77 ± 0.46 × 10⁻² and
21.3 ± 3.3 × 10⁻²), which is the noise of four snapshots, not a bias.

**Computation.**
- A1. Apply the estimator with the dt = 0.005 re-equilibration to *every* stored condition:
  χ ∈ {0.25 … 1.5} × 3 seeds at N = 4000 (`data_paper`), the N = 1000 sweep, N = 9000/16 000
  (`data_box`), N = 22 000 (`data_validate`). Cost: ~2 min per snapshot at N = 4000; 20 snapshots
  × 18 files ≈ 12 CPU-h, or ~4 h with 200 snapshots at the production timestep plus one
  re-equilibrated calibration point per condition. Output: EPR(χ, N) on a single grid.
- A2. Plot ã(χ) and b̃(χ) separately. Decision rule: if b̃/χã is O(1) and χ-independent the
  cross-term is structural screening (expected); if b̃ → const as χ → 0 there is a genuine linear
  term and the positivity argument says the χ = 0 reference is not stationary (aging gel) — which
  is itself the informative outcome.
- Resolves Q1 without the abandoned 13–26 h run, at *all five* system sizes rather than one.

**Manuscript consequence.** §3.10's two-term paragraph and Conclusion 5 are rewritten around
EPR = χ²A(χ) with A rising through unjamming; the "negative cross-term resolved at 5.6σ"
becomes "contact-force screening of ~90% of the drive, resolved at every N". §2.4's caveat about
the subtracted artifact is retired for the excess (it survives only for the raw heat).

## Track B — Frequency-resolved effective temperature via Harada–Sasa (Q3; new cheap data)

**Reasoning result.** The window-dependent T_eff (5.2 → 9.5 → 17.7 kT across T = 200–800) is the
low-frequency divergence of the fluctuation–dissipation violation, not the absence of a
temperature. For overdamped Langevin dynamics with additive noise the Harada–Sasa equality holds
exactly per degree of freedom:

    Q̇_i = (1/s_i)⟨v_i⟩² + (1/s_i) ∫ dω/2π [ C̃_i(ω) − 2T R̃′_i(ω) ]

so the integrated FDT violation over all frequencies *must equal* the dissipation from track A.
That is a closed consistency loop between Q1 and Q3, and T_eff(ω) = C̃(ω)/2R̃′(ω) is the
well-defined Cugliandolo–Kurchan object that §3.10 lacks.

**Why the species-coordinate artifact disappears.** Probe single-particle coordinates, not
collective ones: perturb all particles with a random-sign force f ε_i (ε_i = ±1 i.i.d.) using the
shared-noise protocol already in `fdt.py`; the response of Σ ε_i x_i is N × the mean single-particle
response, cross terms average out, and the coordinates are genuinely independent.

**New data.** From stored equilibrated snapshots at χ ∈ {0, 0.5, 1, 1.5}, N = 4000, 3 seeds:
perturbed and unperturbed twins, t = 2000, positions of 500 tagged particles stored every 0.25
(≈ 64 MB per run). 24 runs × ~2 min. Also gives the velocity autocorrelation at fine resolution,
which nothing in `data_*` currently has.

**Decision rule.** If T_eff(ω) has a plateau over the contact-relaxation band (ω ~ k μ ~ 1–10)
and diverges only below ω ~ 1/t_unjam, tier-II holds at fixed structure and fails only through
structural rearrangement — the boundary §4.7 asserts, now measured. If no plateau exists at any
frequency, Conclusion 8's negative stands and becomes frequency-resolved. Either way the
Harada–Sasa sum rule against track A is a check no other diagnostic in the paper has.

## Track C — Nucleus or crossover (Q4; existing dense data + cheap reruns)

The three `cont2_dense500` files (Δt = 50, 501 snapshots, N = 4000, χ = 1.5) are exactly the data
behind the split/merge table and are enough for C1–C2.

- C1. **Error bars on S\*.** Bootstrap over clusters within seed and over seeds; report the crossing
  with a confidence interval instead of "between the size-5 and size-10 bins". Zero new compute.
- C2. **Size drift ⟨ΔS|S⟩.** Track cluster identity by maximum particle overlap between consecutive
  snapshots. Absorption by the condensate must be treated as a separate hazard h_abs(S), not folded
  into ΔS (it makes the drift look like +3000 at every S). Report the conditional drift for
  non-absorbed clusters and h_abs(S) side by side. A nucleus needs the drift to change sign at S\*;
  the rate crossing does not guarantee that because a split produces two clusters whose sizes
  matter, not just an event count.
- C3. **Committor.** From each dense snapshot, launch 50 noise-resampled runs of t = 500 (10⁴ steps,
  ≈ 30 s each at N = 4000) and record, for every fragment present at t = 0, whether it first reaches
  2S\* (grow), dissolves below S/2 (shrink), or is absorbed. One launch scores ~140 fragments at
  once, so 150 launches (~75 CPU-min) give q(S) with ~10⁴ fragment-fates. A nucleus needs
  q(S\*) ≈ ½ with q monotone in S.
- C4. **Becker–Döring detailed balance and the cycle in size space** (also Q7). With the tracked
  events, test k_split(S+1) n(S+1) ≈ k_merge(S) n(S) n(fragment) size by size. At the late-time
  plateau, measure the net flux around the loop condensate → shed fragment (size S₁) → further
  fragmentation → reabsorption (size S₂ < S₁). If the shedding-size distribution differs from the
  absorption-size distribution, that is a **circulating probability current in the cluster-size
  coordinate** — the first coarse observable in which the irreversibility would project, and it is
  measurable from data already on disk. The χ = 0 null for it comes from track D.

**Manuscript consequence.** Conclusion 7 either gains "critical nucleus, drift sign change at
S\* = 7.4 ± 0.6, q(S\*) = 0.5 ± 0.05" or loses the nucleation analogy for good. C4 may convert
Conclusion 5's "no projection examined" into a positive.

## Track D — A rate-balance model of structure selection (Q5; new cheap data)

The paper has the mechanism qualitatively (activity completes arrested phase separation at small χ,
drives fission at large χ) and the coherence measurement, but no rates as a function of χ and
therefore no prediction of where the dip sits. The dense-snapshot data that make rates measurable
exist only at χ = 1.5.

**New data.** Dense reruns (Δt = 50, t = 2.5 × 10⁴, 500 snapshots) from the stored late
configurations at χ ∈ {0, 0.25, 0.5, 0.75, 1, 1.5}, N = 4000, 3 seeds: 18 runs × ~25 min ≈
7.5 CPU-h, ~1.5 h wall on six cores. These same files serve C4 (χ = 0 null), G1, and F3.

**Analysis.** Per-cluster split rate k_s(S; χ), merge rate k_m(S; χ), cluster propulsion speed
v(S; χ). The minimal model is coagulation–fragmentation with k_m ∝ v(χ) ∝ χ (measured, not
assumed: §3.5 says clusters move 10.7× faster at matched size) and k_s(χ) rising from zero above a
threshold where internal nonreciprocal stress exceeds spring cohesion. Its stationary n_cl(χ) has a
minimum where the two rates cross. **Decision rule:** the model predicts the dip at χ_min from
rates measured at χ ≠ 0.25 alone; the measured dip (2.2 clusters at χ = 0.25) is the test.

**Manuscript consequence.** §4.7's "no functional form predicts the dip" becomes "structure is a
rate-balance fixed point, and here is the form". This is also the sharpest answer available to
NWAP's structural question: the selected structure is a dynamical fixed point of event rates, not
an extremum of a static functional, which is exactly the boundary §4.7 draws by assertion.

## Track E — Coarsening law from the structure factor (Q6; existing data, zero compute)

Every stored run has 100–500 snapshots of positions. The standard active-matter analysis —
domain size L(t) from the first moment of S(k), or equivalently the inverse of the pair-correlation
crossing — has not been done; the paper reasons from cluster counts and lcf alone.

- E1. L(t) at every χ and N. Growth exponent 1/z vs χ. Arrest in the conventional sense is
  logarithmic or saturating L(t); transience is a power law that only the box cuts off.
- E2. Finite-size collapse of lcf(t) and L(t) across N = 4000, 9000, 16 000, 22 000 against
  t/N^{z/2}. This *measures* z instead of extrapolating 1/z ∈ [2.4, 3.7] from the fourfold range
  of N, and gives the transient lifetime scaling directly.
- E3. At χ = 0 the same curve is the gel's aging law, and its comparison with χ = 0.25 quantifies
  unjamming as a change of growth exponent rather than as a cluster-count ratio.

**Decision rule.** Power-law L(t) with 1/z consistent across N ⇒ Conclusion 7's conditional
becomes measured; saturation that does not move with N ⇒ a genuine selected scale, and §3.8 is
wrong. Either result is a headline for the transience claim.

## Track F — The directed contact graph and what circulation is allowed (Q7, Q8; existing data)

- F1. **Build the directed graph.** Because `lvec` is bi-valued, the pair propulsion on every L–S
  contact points from S to L and vanishes on L–L and S–S contacts; the edge flow is w_ij = (a_j−a_i)/2,
  antisymmetric, with magnitude depending on r_ij. Hodge-decompose it with the `hodge_dominance.py`
  machinery. Prediction from the algebra: the flow is a two-level potential (φ_L − φ_S) up to the
  r-dependence of the kernel, so the gradient fraction should exceed ~0.9 and the curl residual
  should scale with the variance of the kernel over the contact shell. Measure it and its null
  (weights shuffled among edges). This closes the Limitation exactly as stated, with the caveat
  §4.2 already makes (edge-flow Hodge ≠ configuration-space current) carried over verbatim.
- F2. **Parity.** The model is achiral: the force law, noise and geometry are invariant under
  reflection. Every spatial pseudoscalar — cluster angular velocity, vorticity of the coarse
  velocity field, cluster spin — therefore has zero ensemble mean at every χ *by symmetry*, not by
  measurement. The observable-plane signed-area rates of §2.4 are *not* pseudoscalars under
  reflection and are not covered; they remain empirical nulls. This one-paragraph argument belongs
  in §4.2 and explains why the chiral Vicsek control (§3.11) circulates while this system cannot:
  the property being measured is a broken parity, not a broken time reversal.
- F3. **All-pairs lag asymmetry.** On the full set of scalar time series (n_cl, lcf, turnover,
  σ²_v, U, static EPR from A1), compute the antisymmetric part of the lagged covariance matrix
  C(τ) − C(−τ) at every lag against the χ = 0 null. This is the Gladrow–Broedersz test in its
  strongest form on the observables already computed; if it is null too, "no projection examined"
  is upgraded to "no cyclic mode among these observables at any lag".

## Track G — Cheap closures of stated limitations

- G1. σ²_v as a function of sampling interval from the dense runs of track D, plus one dense
  monodisperse rerun; report the interval at which the mono/bi ratio reaches the source's 4. Closes
  the "cannot be confirmed" limitation to the extent possible without the source interval.
- G2. iAAFT surrogates for the cilia analysis (data on disk, `cilia_analysis.py` takes a surrogate
  function). Report z-scores under both nulls. If the paper's numbers survive, say so; if they
  halve, say that.
- G3. Green–Kubo at χ = 0 (Q9). With track B's fine-resolution χ = 0 twins, compute
  L_aa = (1/2T)∫⟨J_a(t)J_a(0)⟩ dt with J_a = F_a∘ẋ; the t = 0 delta term is track A's ã, and the
  remainder is the relaxation correction. L_aa is the linear-response prediction of
  lim_{χ→0} EPR/χ². Compare with A1's measured A(χ) at χ = 0.25. This is the "defensible test" §3.10
  says it lacks; it is defensible because J_a is time-odd, so Onsager symmetry of L is a theorem here
  and the only question is whether the measured low-χ curvature matches the equilibrium
  correlation. The χ = 0 reference is an aging gel, so agreement is expected only for t < t_unjam(χ),
  which track E measures.
- G4. Rectangular box at N = 22 000: ~5 CPU-days at the source duration. Not worth it unless a
  referee asks; the 3% at N = 4000 stands.

## Track H — Where new physics or new data are genuinely required

- H1. **Hydrodynamic interactions.** §5 names them as "exactly the class of ingredient that could
  stabilise finite clusters". Testing that means a different model. The minimal honest version: add
  a reciprocal far-field mobility coupling (a 2D near-wall Stokeslet, decaying as 1/r² with the
  source geometry's wall) and repeat track E at N = 4000 for χ = 1.5. If L(t) saturates
  N-independently with hydrodynamics and not without, the source experiment's arrest is
  hydrodynamic and the agent-based model cannot show it. ~2 CPU-days once the kernel exists; a
  separate short paper if it works.
- H2. **Thermodynamic-limit lifetime.** Track E measures z from four sizes; a fifth at N = 40 000
  (L = 76) to twice the current duration would be the check. ~3 CPU-days per seed on the M5 Studio.
- H3. **Dominance power.** Only raw interaction tracking resolves the cyclicity question;
  `REQUEST_RAW_TRACKING.md` already specifies it. No computation can substitute.
- H4. **Source sampling interval for σ²_v.** One email to the authors of ref. [1]; G1 bounds it
  meanwhile.

## Execution order and what each step feeds

1. A1 with dt calibration (running now) → single-N and all-N EPR grid. Feeds B (sum rule), G3.
2. Track D dense reruns launched in parallel (1.5 h wall) → feeds C4, D, F3, G1.
3. C1, C2, E1–E3, F1, F3, G2 on existing data while D runs: no compute contention.
4. C3 committor launches (75 CPU-min) after D finishes.
5. B twins (24 × 2 min), then G3.
6. Write-up: §3.10 rebuilt around A + B, §3.9 around C, §3.5/§4.7 around D, §3.8 around E,
   §4.2/§5 around F.

Everything above through step 5 fits in one working day of wall time on the M2 Max, and none of it
depends on the outcome of any other step except as noted.
