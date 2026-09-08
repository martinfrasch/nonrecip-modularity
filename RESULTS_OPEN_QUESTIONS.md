# Results log for OPEN_QUESTIONS_PLAN.md

Running record, newest at the bottom of each track. Numbers here are the ones to carry into the
manuscript once each track closes; anything marked *provisional* is awaiting a calibration or
a replicate.

## Track A — dissipation from configurations (Q1, Q2)

**Estimator.** `static_epr.py`. Per-particle excess dissipation EPR/N = χ(b̃ + χã) with
ã = Σ|F_a,i|²/(s_i T N) and b̃ = [Σ F_a,i·F_s,i/s_i + Σ D_i ∇_i·F_a,i]/(T N), both averaged over the
late half of every stored run at the production timestep. `epr_static.csv` holds all 107 files.

**Validation against the paired trajectory probe** (N = 4000, three seeds, 100 late snapshots each):

| χ | paired probe (paper) | configurational | ratio |
|---|---:|---:|---:|
| 1.5 | 6.53e−3 (N = 1000) | 6.26 ± 0.42e−3 | 0.96 |
| 2 | 1.22e−2 | 1.22 ± 0.10e−2 | 1.00 |
| 3 | 2.63e−2 | 2.34 ± 0.16e−2 | 0.89 |
| 5 | 7.25e−2 | 6.46 ± 0.34e−2 | 0.89 |
| 8 | 1.88e−1 | 1.78 ± 0.10e−1 | 0.94 |

**Timestep bias.** The symmetric-sector identity Q̇_s = −dU/dt fails at dt = 0.05 (Q̇_s ≈ +4.4 T
per particle per unit time against a measured −dU/dt ~ 10⁻⁶) and recovers under re-equilibration
at dt = 0.005 (Q̇_s ≈ 0.2–0.6). The calibrated excess at χ = 2 is 1.52e−2 on four snapshots
against 1.24e−2 uncalibrated, i.e. the production-timestep values may sit ~10–20 % low at high
drive — the same direction and size as the ratios above. *Provisional: a trajectory-mode
calibration (100 small-dt configurations per file) is queued for the low-drive files, where the
question is whether the excess is zero.*

**The single-N low-drive grid** (N = 4000, L = 24, base runs and cont1 only; the abandoned
13–26 h probe was to produce exactly this row):

| χ | 0.25 | 0.5 | 0.75 | 1 | 1.5 | 2 | 3 | 5 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| EPR/N (×10⁻³) | 0.02 ± 0.13 | 0.05 ± 0.05 | 0.58 ± 0.06 | 2.26 ± 0.12 | 6.26 ± 0.42 | 12.2 ± 1.0 | 23.4 ± 1.6 | 64.6 ± 3.4 | 178 ± 10 |
| EPR/(Nχ²) (×10⁻³) | 0.3 | 0.2 | 1.0 | 2.3 | 2.8 | 3.1 | 2.6 | 2.6 | 2.8 |
| screening −b̃/(χã) | 0.99 ± 0.13 | 0.99 ± 0.02 | 0.96 ± 0.01 | 0.91 ± 0.01 | 0.89 ± 0.02 | 0.88 ± 0.02 | 0.90 ± 0.01 | 0.91 ± 0.01 | 0.91 ± 0.01 |

**Reading.**
1. ã, the naive |F_a|² dissipation, is nearly constant in χ (0.025–0.030): it is set by the
   number of L–S contacts, not by the drive.
2. b̃ is proportional to χ at every drive with coefficient −0.89 ã above χ ≈ 1. The
   "negative linear term" of §3.10 is therefore a *quadratic* structural screening term: about
   nine tenths of the nonreciprocal push is balanced by contact forces and does not dissipate.
   As reasoned in the plan, a genuine linear term with b < 0 is forbidden by EPR ≥ 0.
3. Below χ ≈ 0.75 the screening is complete to within error and the excess dissipation is
   consistent with zero. EPR/χ² is not a constant with a linear correction; it is a **sigmoidal
   onset** between χ = 0.5 and χ = 1.5 that plateaus at 2.7 × 10⁻³. The onset coincides with
   the structural transition of §3.5 (n_cl minimum at χ = 0.25, fragmentation above): in the
   condensed, jammed state the drive is stored as contact stress; dissipation requires relative
   motion, which begins when the antisymmetric push exceeds what the contacts can hold.
4. At χ = 0 the identity b̃ = 0 holds to within noise (mean +0.0003 on ã = 0.030), so the
   kinetically arrested gel is locally Boltzmann along the F_a direction; the χ > 0 screening is a
   drive-induced structural response and not a property of the reference.
5. Across system size at χ = 1.5: EPR/(Nχ²) = 3.7 (N = 1000), 2.8 (4000), 2.4 (9000),
   2.4 (16 000) × 10⁻³. The N = 1000 pilot value is ~30 % high, which is the direction in which
   the N = 1000 two-term fit of §3.10 was pulled.

## Track C — nucleus or crossover (Q4)

`cluster_kinetics.py` on the three `cont2_dense500` files (χ = 1.5, N = 4000, Δt = 50):
209 748 cluster-observations of size ≥ 2 excluding the condensate.

**Event definitions differ from the manuscript's (unsaved) analysis** and the difference matters:
fission (≥ 2 pieces of size ≥ 2) is separated from monomer evaporation, and fusion with another
size-≥ 2 cluster from monomer capture. With that separation the fission probability rises
monotonically with size (0.03 at S ≈ 3, 0.24 at S ≈ 8, 0.34 at S ≈ 23, 0.54 above 100) and
fusion is flat at 0.17–0.26, so **there is no split/merge crossing** in this bookkeeping; the
manuscript's crossing at 7–8 came from counting evaporation as splitting.

**The Becker–Döring drift changes sign.** Restricting to monomer/dimer exchange events (no
fission, no fusion), the size drift per Δt = 50 is

| S bin | 2 | 3–4 | 5–6 | 7–9 | 10–13 | 14–19 | 20–29 | 30–49 | 50–99 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ⟨ΔS⟩ | +0.092 | −0.139 | −0.161 | −0.056 | +0.022 | +0.045 | +0.092 | +0.123 | +0.143 |
| sem | 0.003 | 0.004 | 0.009 | 0.012 | 0.017 | 0.022 | 0.032 | 0.041 | 0.055 |

negative from S = 3 to S ≈ 9 and positive above. **Drift zero at S\* = 10.1, 68 % block-bootstrap
interval 9.7–10.7, defined in 100 % of resamples.** That is the unstable fixed point the
manuscript said it had not measured. Dimers are the exception (they gain on net). The
committor (C3) remains to be run; the plan's decision rule is q(S\*) ≈ ½.

**Absorption as a hazard** is small and flat: 1.5–3.5 % per interval at every size. Including
it in the drift gives +55 to +140 at every S (condensate size), which is the artifact
FOLLOWUP.md T2 warned about.

**No cycle in size space (C4, χ = 1.5).** Over the late plateau the condensate shed 4 435 pieces
(mean size 5.60, median 3) and absorbed 5 114 clusters (mean 5.49, median 3); a two-sample KS test
on the two size distributions gives D = 0.010, p = 0.97. What is shed and what is reabsorbed are
the same population, so the shedding–reabsorption loop carries no detectable circulating current
in the cluster-size coordinate. Conclusion 5's "no projection examined" survives in this
coordinate too. The χ = 0 comparison waits on the track D reruns.

## Track D — rate balance across χ

Dense reruns (Δt = 50, t = 2.5 × 10⁴) for χ ∈ {0, 0.25, 0.5, 0.75, 1} × 3 seeds plus one
monodisperse seed are running in `data_dense/` (16 jobs, five workers). Analysis follows with
the same script.

## Track E — coarsening law

`coarsening_law.py` running over all chained runs; output `coarsening.csv`.

## Track F — directed contact graph (Q8) and what circulation is allowed (Q7)

**F1, Hodge decomposition of the antisymmetric edge flow** (`hodge_contact.py`, six late snapshots
per file, N = 4000). The flow w_ij = (a_j − a_i)/2 on contact edges, decomposed on the contact
graph itself into a node potential and a cycle-space residual:

| χ | 0 | 0.25 | 1 | 1.5 | 8 |
|---|---:|---:|---:|---:|---:|
| cyclic fraction ‖w − ∇s‖/‖w‖ | 0.077 | 0.079 | 0.075 | 0.071 | 0.063 |
| permuted-magnitude null | 0.738 | 0.742 | 0.729 | 0.735 | 0.740 |
| L–S-support null | 0.749 | 0.747 | 0.739 | 0.741 | 0.749 |
| CV of \|w\| on L–S edges | 0.12 | 0.13 | 0.12 | 0.12 | 0.11 |

The antisymmetric sector of this model, as an edge function on the contact graph, is **99.4 %
gradient in squared norm at every χ**, against 45 % for a random flow on the same graph. The
residual 6–8 % is what the 11–13 % spread of the kernel over the contact shell predicts: the flow
is a two-level potential (S → L on every mixed contact) up to the r-dependence of g. This is the
directed-graph test the Limitations section named, and it comes out on the side of Conclusion 11:
the antisymmetric coupling of this model is not circulating, it is a gradient, and the reason no
coarse observable circulates is that nothing in the pair law does. The caveat of §4.2 stands: this
is a Hodge decomposition of a static edge function, not of the configuration-space current.

**F2, parity.** The force law, noise and geometry are reflection-invariant, so every spatial
pseudoscalar (cluster spin, vorticity of the coarse velocity field) has zero mean at every χ by
symmetry. The observable-plane signed-area rates of §2.4 are not pseudoscalars and are not
covered; they remain empirical nulls. One paragraph for §4.2.

**F3** waits on the dense reruns (the paper-scale series has 100 late points per seed).

## Track G — cheap closures

**G2, iAAFT surrogates for the cilia comparison** (`cilia_analysis.py --iaaft`, 30 surrogates
per axoneme, 184 axonemes). The Limitations section predicted the stricter null "would likely
reduce the reported z-scores". It barely does:

| null | median \|z\| | frac \|z\| > 2 | frac \|z\| > 3 |
|---|---:|---:|---:|
| phase randomisation (paper) | 3.2 | 0.92 | 0.55 |
| iAAFT (spectrum + amplitude distribution) | 3.0 | 0.92 | 0.51 |

Per-ATP rows move by at most 0.3 in median |z|. The cilia circulation result stands under the
amplitude-preserving null; the limitation can be closed with the second row.

### Track E results (Q6) — `coarsening.csv`, `coarsening_fits.csv`

Chained base + continuation runs, every second snapshot, three observables per snapshot: the
first-moment length L_k of S(k) (k ≤ π), the mass-weighted mean cluster size S_w = ΣS²/ΣS, and lcf.

1. **The condensation time scales linearly with N.** Time to lcf = 0.5 at χ = 1.5: 5.5 × 10⁴
   (N = 4000), 1.44 × 10⁵ (9000), 2.52 × 10⁵ (16 000), i.e. t_½ ∝ N^{1.0–1.2}; at χ = 1: 5.7 × 10⁴
   and 1.26 × 10⁵ for N = 4000 and 9000 (N^0.97). The lcf(t) curves collapse in t/N and not in
   t/√N:

   | t/N | 5 | 10 | 20 | 40 |
   |---|---:|---:|---:|---:|
   | lcf, N = 4000 | 0.30 | 0.47 | 0.55 | 0.77 |
   | lcf, N = 9000 | 0.30 | 0.44 | 0.60 | 0.77 |
   | lcf, N = 16 000 | 0.21 | 0.32 | 0.55 | 0.75 |

   So the transient's lifetime is t_x ∝ L², a diffusive crossing of the box, *measured* over a
   fourfold range of N. Conclusion 7's conditional "N^(1/z) with 1/z between 2.4 and 3.7" should
   be replaced by this.
2. **Growth law.** In the growth phase (t = 10⁴–2 × 10⁵) S_w grows as t^{0.5–1.0} at χ = 1–1.5 and
   saturates (exponent ≈ 0.05) once the condensate has formed; the largest N = 22 000 replicate
   shows S_w ∝ t^{0.55} over 95 points with no saturation before lcf = 0.5. A power law that only
   the box cuts off is a transient, not an arrested state — the conventional test, and it agrees
   with §3.8.
3. **Unjamming is a change of growth exponent.** At χ = 0 the gel ages with S_w ∝ t^{0.11} and
   never reaches lcf = 0.5; at χ = 0.25 the same initial ensemble grows with exponent 0.58 and
   condenses by t = 4 × 10⁴. Above χ = 3 the exponent falls again (0.43 at χ = 5, 0.20 at χ = 8,
   fitted on ≤ 9 pre-saturation points, so provisional): strong drive slows coarsening in the
   growth phase even though the condensate still forms. That is the "arrested coarsening" of
   the source paper appearing as a *rate*, not as a selected size.

### Track D results (Q5) — dense reruns at every χ, `kinetics_events.csv`

Δt = 50, t = 2.5 × 10⁴, three seeds per χ from the stored late configurations; the χ = 1.5 row is
the existing `cont2_dense500` set. Fragments = clusters of size ≥ 2 other than the condensate.

| χ | 0 | 0.25 | 0.5 | 0.75 | 1 | 1.5 |
|---|---:|---:|---:|---:|---:|---:|
| fragments per snapshot, observed | 9.2 | **1.1** | 15.0 | 57 | 91 | 140 |
| condensate shedding rate (pieces per Δt) | 0.063 | 0.088 | 0.70 | 2.15 | 3.23 | 2.96 |
| reabsorption probability per fragment per Δt | 0.0073 | **0.082** | 0.047 | 0.039 | 0.038 | 0.023 |
| **fragments predicted = shedding / reabsorption** | 8.6 | 1.1 | 14.9 | 55 | 86 | 129 |
| fission probability, S = 10–29 | — | — | 0.10 | 0.21 | 0.25 | 0.30 |
| Becker–Döring drift zero S\* | — | — | 7.3 | 10.7 | 9.3 | 10.1 |
| speed of clusters S ≥ 30 (×10⁻⁴) | 0.9 | 0.8 | 2.0 | 2.2 | 2.9 | 3.7 |
| shed vs absorbed size, KS p | 1.00 | 1.00 | 0.99 | 0.86 | 0.53 | 0.97 |

**The functional form that predicts the dip.** At every χ, including the reciprocal reference,
the fragment count is the condensate's shedding rate divided by the per-fragment reabsorption
probability, to within 8 %. That is a two-rate balance, n_frag = k_shed / p_abs, and the
non-monotonic curve of §3.5 is the product of two monotone but *differently thresholded* rates:
reabsorption switches on at χ = 0.25 (×11 over χ = 0, because weak activity mobilises fragments
and the condensate) while shedding is still at its reciprocal value, so fragments are cleared
faster than they are made and n_cl falls fivefold; shedding then switches on between χ = 0.25 and
0.5 (×8) and again to 0.75 (×3), and n_frag rises by two orders of magnitude. The fission
probability at fixed size follows the same onset (zero at χ ≤ 0.25, 0.10 at 0.5, 0.30 at 1.5).
Structure is a **rate-balance fixed point**, and no static functional is needed to predict it.

**The shedding onset and the dissipation onset coincide.** Track A puts the onset of excess
dissipation between χ = 0.5 and 1, where the contact screening drops from 0.99 to 0.91; the
fission and shedding onsets sit between 0.25 and 0.75. Both are the same threshold: the
antisymmetric push exceeding what the contact network can hold statically. Below it the drive is
stored as stress and neither fragments nor dissipates; above it it does both.

**S\* is χ-independent above the onset**: 9–11 particles at χ = 0.75–1.5 (7 at χ = 0.5, with a
tenth of the events). It is a property of the contact energetics, not of the drive, as §3.9
guessed from the N- and density-independence.

**C4 at every χ.** The shed and reabsorbed size distributions are indistinguishable at all six
χ (KS p ≥ 0.53). There is no circulating current in the cluster-size coordinate at any drive.

### Track A, calibrated (`epr_static_calibtraj.csv`)

Trajectory-mode calibration: from two late snapshots per file, 500 time units at dt = 0.005,
100 configurations each, first fifth discarded (162 configurations per file, three seeds per
χ). The stationarity residual Q̇_s falls from 4.4 to 0.35. **These supersede the production-
timestep row above at low drive**, where the excess is a small residual of two cancelling terms
and a 5 % shift in the screening term moves it by a factor of three.

| χ | 0.25 | 0.5 | 0.75 | 1 |
|---|---:|---:|---:|---:|
| EPR/N, production dt (×10⁻³) | 0.02 | 0.05 | 0.58 | 2.26 |
| **EPR/N, calibrated (×10⁻³)** | **−0.03 ± 0.02** | **0.24 ± 0.03** | **1.61 ± 0.17** | **2.32 ± 0.22** |
| EPR/(Nχ²), calibrated (×10⁻³) | −0.5 | 0.95 | 2.87 | 2.32 |
| screening −b̃/(χã), calibrated | 1.02 | 0.97 | 0.90 | 0.91 |

The onset is sharper and earlier than the uncalibrated row suggested: zero at χ = 0.25 (three
seeds, each within 1σ of zero), a third of the plateau at χ = 0.5, and the plateau of
2.3–2.9 × 10⁻³ reached by χ = 0.75. The screening fraction steps from 1.0 to 0.90 across the same
interval and is flat at 0.89–0.91 from χ = 0.75 to χ = 8. This is the same threshold as the
shedding and fission onset of track D (0.25 → 0.5 → 0.75). The two-term fit of §3.10 was
describing this step with a negative linear term; the physical description is a threshold in χ
below which the contact network holds the antisymmetric push statically.

## Track B — effective temperature on single-particle coordinates (Q3)

`hs_fdt.py` + `hs_analysis.py`: random-sign force f ε_i x̂ on every particle (f = 10⁻⁵, linear to
t ≥ 200 against f = 10⁻⁶), shared noise stream with an unperturbed twin, T = 1000, positions every
1.0. Twelve twins per χ (three seeds × four late starts), N = 4000. The diagonal response
χ(t) = ⟨ε_i Δx_i⟩/f is conjugate to the single-particle MSD, so T_eff(t) = MSD/2χ is a proper
Einstein ratio per degree of freedom — the object the species-coordinate measurement of §3.10
could not provide. `hs_results.csv`.

| χ | T_eff/T at t = 1–10 | t at which T_eff/T exceeds 1.1 / 1.2 | T_eff/T at t = 100 / 400 / 800 | L-species ÷ S-species at t = 400 / 800 |
|---|---:|---:|---:|---:|
| 0 | 1.00–1.04 | never / never | 1.00 / 0.98 / 0.99 | 1.00 / 1.02 |
| 0.5 | 1.02–1.04 | 216 / 340 | 1.06 / 1.28 / 1.74 | 1.01 / 0.98 |
| 1 | 1.01–1.05 | 37 / 74 | 1.33 / 2.60 / 3.94 | 1.03 / 0.93 |
| 1.5 | 1.05–1.08 | 13 / 30 | 1.81 / 4.38 / 7.17 | 1.02 / 0.89 |

Three findings, in order of strength.

1. **The estimator calibrates exactly.** At χ = 0 the ratio is unity within 3 % at every lag from
   1 to 800 on both species (twelve twins). The 0.21 and 0.11 that §2.5/fdt.py reported for
   single-particle coordinates were the non-conjugate pairing of a species-wide force with a
   single-particle MSD, not a property of the system.
2. **There is a single effective temperature, and it equals T, over a drive-dependent window.**
   FDT holds to within 10 % up to t ≈ 220 at χ = 0.5, 37 at χ = 1 and 13 at χ = 1.5, i.e. for
   ω ≳ 1/t_FDT(χ). Below that frequency the ratio grows as (T_eff − 1) ∝ t^{1.0} at χ = 1 and 1.5
   (the paper's t^0.89 on the species coordinate is this regime), without a plateau to t = 800.
3. **The two species share the same T_eff(t) at every lag**, to 2 % at t = 400 and 10 % at
   t = 800, despite mobilities differing by 1.5×. The single-temperature question §3.10 said it
   could not ask — because its two coordinates were one degree of freedom — has an answer on
   independent coordinates: one temperature, shared, but a function of the observation time,
   not a number. Whether that is "tier II at fixed structure" is a matter of definition; what is
   measured is FDT with T_eff = T above a cutoff frequency that falls as 1/χ² (t_FDT: 216, 37, 13
   for χ = 0.5, 1, 1.5), and a shared, growing violation below it.

The Harada–Sasa sum rule was not evaluated: a frequency-domain estimate from second differences
of the MSD is too noisy at this sampling, and the time-domain result above carries the
conclusion. The violation lives entirely below ω ≈ 1/t_FDT(χ), which places the dissipation of
track A in cluster-scale motion rather than at contact scale.
