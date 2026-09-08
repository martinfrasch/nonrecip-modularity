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
