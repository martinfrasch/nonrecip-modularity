# Is a finite characteristic cluster size selected?

Status as of 2026-08-26. **Two of three points converged and agreeing; the decisive
composition-corrected replicate is still running.** Treat the conclusion as provisional.

## The question

Claim 1 of the NWAP reading (`NWAP_ASSESSMENT.md`) requires the system to *select* a finite
cluster size — an interior peak in the cluster-size distribution. At our original box (N=4000,
L=24) the steady state is instead one system-scale condensate holding ~78% of particles plus a
monotonically decaying spray of dimers and trimers, with no interior peak. But a condensate that
spans the box is what a too-small system shows regardless, so this could not be decided there.

**The test:** hold density fixed and grow the box. A finite characteristic size S* forces the
largest cluster to stay constant (lcf ∝ S*/N, falling as 1/L²). True phase separation makes the
largest cluster grow ∝ N (lcf constant). The exponent in `largest ~ N^x` distinguishes them:
x=0 for a characteristic size, x=1 for phase separation.

Densities are matched exactly at 6.944 particles per unit area: (N=4000, L=24), (N=9000, L=36),
(N=16000, L=48).

## Convergence is the dominant difficulty

Larger boxes coarsen far more slowly, and at fixed t̂ the larger boxes are systematically less
equilibrated — which *depresses* the largest cluster and mimics a characteristic size. Every
continuation run so far has moved the largest cluster **up**:

| condition | window 1 (t̂=0–4e5) | window 2 (t̂=4e5–8e5) | change |
|---|---:|---:|---:|
| N=4000, χ=1.5 | lcf 0.753 ± 0.040 | lcf 0.779 ± 0.004 | +3.4% — converged |
| N=9000, χ=1.5 | lcf 0.737 ± 0.129 | lcf 0.802 ± 0.001 | +8.8% — converged |
| N=9000, χ=1 | lcf 0.443 ± 0.005 | lcf 0.827 ± 0.006 | **+87% — was nowhere near steady state** |

The χ=1 case at L=36 is the cautionary one: at t̂=4e5 it read as a capped cluster size, and
doubling the run nearly doubled the condensate. Any conclusion drawn from an unconverged large
box will falsely favour a characteristic size. Note also that seed scatter collapses on
convergence (±0.129 → ±0.001), which is a more reliable convergence indicator here than a
slope fit — cluster counts fluctuate ~19% within a run in the crossover region.

## Current result (χ=1.5)

| N | L | largest cluster | ratio | (N ratio) | lcf | status |
|---:|---:|---:|---:|---:|---:|---|
| 4,000 | 24 | 3,116 | 1.00× | 1.00× | 0.779 | converged (t̂=8e5) |
| 9,000 | 36 | 7,221 | 2.32× | 2.25× | 0.802 | converged (t̂=8e5) |
| 16,000 | 48 | 8,335 | 2.68× | 4.00× | 0.521 | **not converged** (t̂=4e5) |

- **Exponent from the two converged points: 1.04** — the phase-separation value, against 0.00
  for a finite characteristic size.
- lcf is flat to slightly rising as the system more than doubles (0.779 → 0.802).
- Three-point fit reads 0.73, but that is produced entirely by the unconverged N=16000 point.
- Fragment population is extensive: n_cl ∝ N^0.96, median cluster size 3 at every box.

**Reading: no finite characteristic cluster size, at this composition.** The condensate grows in
proportion to the system, with a system-proportional spray of small fragments alongside it.

## The caveat that outranks all of the above

Every run in the table uses **25% large particles**. The source paper uses **5000/22000 =
22.7%** (`NWAP_ASSESSMENT.md`, validation section). Large particles carry 5.1× the EHD strength
(l⁴ scaling), so an excess of them biases toward condensation — precisely the effect being
measured. Our packing fraction is also low (35.4% vs the paper's 38.0%).

A replicate at the paper's own specification is running: N=22,000, 22.7% large, t̂=3.6e5
(=3600 s, their duration), equal-area square box L=53.67. Until it lands, the conclusion above
describes *our parameterisation*, not necessarily the published system.

The paper reports that in the head-large geometry a significant fraction of particles remain in
small clusters, whereas our runs put ~78% in one condensate. Our equal-size and tail-large cases
*do* reproduce the paper's contrast (lcf 0.94–1.00, a single cluster), so the qualitative
contrast holds; it is the head-large mass distribution that is unverified.

## Deviations from the paper in the validation replicate

- **Square box, not 1.8:1.** The paper's domain is 648×360 μm (72×40 λ). Our kernel assumes a
  square box, so the replicate uses an equal-area square (L=53.67). Cluster sizes are small
  compared to either dimension, so this should be second-order, but it is a deviation and not
  a fix.
- **σ²_v is not comparable.** The paper reports mono ≈ 5e-3, bidisperse ≈ 2e-2 μm²/s² (ratio 4×);
  ours gives a ratio of ~90×. Ours is the variance of *speed* coarse-grained over one snapshot
  interval (10 s), not an instantaneous velocity variance — see `AUDIT.md`. The discrepancy is
  almost certainly this definitional difference rather than a model mismatch, but it cannot be
  confirmed without the paper's sampling interval.

## Reproduce

```bash
python simulate.py --sweep chibox --chi-list 1.5 --seeds 2 --N 9000  --box 36 --T 4e5 --nsnap 200
python simulate.py --sweep chibox --chi-list 1.5 --seeds 2 --N 16000 --box 48 --T 4e5 --nsnap 200
python continue_run.py --pattern "data_box/*N9000*.npz" --T 4e5      # convergence is essential
python simulate.py --sweep validate --seeds 2 --N 22000 --box 53.666 \
       --frac-large 0.227273 --T 3.6e5 --nsnap 200                   # paper specification
```
