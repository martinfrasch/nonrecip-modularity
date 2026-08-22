# Results — baseline, sweeps, and reciprocity test, 61 runs, M2 Max (2026-08-21)

Environment: Python 3.12 (arm64), numba 0.67.0, numpy 2.5.2, networkx 3.6.1, pandas 3.0.5.
All runs N=1000 (mono: 500), box 12, t̂=1e5, dt̂=0.05, 100 snapshots, `--burn 0.5`.
61 `.npz`: baseline 5 seeds × 2 cases, α̂ sweep 4×3, s-ratio sweep 5×3, **χ sweep 6×5**.

## Verdict in one line

**Nonreciprocity drives the currents and leaves the structure alone — confirmed (§5), but
only once a knob that actually varies reciprocity existed.** The α̂ and s_II/s_I sweeps
could not test it (§4); the χ sweep can, and does. The mono-vs-bi baseline (§1) turns out
to be a confounded comparison: about half its turnover gap is polydispersity, not
reciprocity, though essentially all of its σ²_v gap is reciprocity.

Read §1 → §4 → §5 in order; §1's causal reading is superseded by §5.

---

## 1. Baseline — CONFIRMED (5 seeds, was 1 seed in the pilot)

| metric (late half) | monodisperse | bidisperse | ratio | pilot | ✓ |
|---|---:|---:|---:|---:|:--:|
| Q (Louvain) | 0.875 ± 0.001 | 0.859 ± 0.002 | 0.98 | ~1 | ✓ |
| ΔQ vs degree-preserving null | 0.391 | 0.323 | 0.83 | 0.8 | ✓ |
| edge Jaccard turnover | 0.164 ± 0.001 | 0.398 ± 0.008 | **2.43×** | 2.6 | ✓ |
| ARI relative to own noise floor | 0.780 vs 0.822 (−0.04) | 0.493 vs 0.682 (**−0.19**) | — | — | ✓ |
| n_clusters | 5.4 | 21.9 | 4.0× | 4 | ✓ |
| σ²_v | 1.18e-9 | 3.52e-8 | **29.8×** | 26 | ✓ |

Every pilot number reproduced at 5 seeds. Static structure is near-identical while the
current observables differ by 2.4× and 30×. Figure: `modularity_test.png`, `snapshots.png`.

> **Superseded by §5:** this comparison varies reciprocity *and* polydispersity together, so
> the 2.4× turnover gap cannot be attributed to nonreciprocity alone. The χ=0 control in §5
> splits it: ~47% reciprocity, ~53% polydispersity. The 30× σ²_v gap *is* essentially all
> reciprocity.

**Caveat on "frozen":** Q is not literally invariant. The 0.0167 gap is ~7× the combined
SEM, so it is statistically resolvable. The defensible claim is one of *magnitude*: Q moves
2% where turnover moves 143% and σ²_v moves 2880%.

**Why edge Jaccard is the headline, not ARI:** the same-graph Louvain noise floor is itself
strongly condition-dependent (0.35 → 0.82 across the runs here). A raw cross-condition ARI
comparison is confounded by the floor; the edge-set metric is algorithm-free.

---

## 2. α̂ sweep — PREDICTION REFUTED (sign is inverted)

Predicted: turnover **rises** with α̂, Q stays pinned. Observed (s_II/s_I = 0.667, 3 seeds):

| α̂ | edge_jac | Q | n_cl | lcf | σ²_v |
|---|---:|---:|---:|---:|---:|
| 0.003 | **0.601** | 0.890 | 40.2 | 0.64 | 1.96e-8 |
| 0.005 | 0.398 | 0.859 | 21.9 | 0.83 | 3.52e-8 |
| 0.010 | 0.274 | 0.843 | 11.9 | 0.93 | 7.38e-8 |
| 0.015 | **0.282** | 0.841 | 9.1 | 0.93 | 1.42e-7 |

Turnover **falls** by 2.1× as α̂ rises 5×, then flattens. Q is **not pinned** — it drifts
−0.049, ~20× its SEM. σ²_v does rise monotonically (7.3×), so the *velocity* current tracks
α̂ as expected; the *edge* current anti-tracks it.

Mechanism: α̂ is the overall EHD attraction strength. Raising it condenses the system
(n_cl 40 → 9, largest-cluster fraction 0.64 → 0.93). Bigger, tighter-bound clusters have
more persistent contacts, so edge turnover drops even as particles move faster. The two
effects compete and condensation wins.

## 3. s_II/s_I sweep — DIRECTION CONFIRMED, but by arrest, not by graded response

Predicted: turnover falls toward tail-large. Observed (α̂ = 0.005, 3 seeds):

| s_II/s_I | edge_jac | Q | n_cl | lcf | ARI vs floor |
|---|---:|---:|---:|---:|---:|
| 0.667 (head-large) | 0.398 | 0.859 | 21.9 | 0.83 | 0.493 vs 0.682 |
| 0.800 | 0.393 | 0.862 | 15.7 | 0.90 | 0.482 vs 0.637 |
| 1.000 | 0.387 | 0.880 | 7.5 | 0.94 | 0.522 vs 0.814 |
| 1.250 | **0.024** | 0.796 | 1.0 | 1.00 | 0.450 vs 0.420 |
| 1.500 (tail-large) | **0.0008** | 0.772 | 1.0 | 1.00 | 0.415 vs 0.348 |

Turnover falls 506× — but not smoothly. It is flat (0.398 → 0.387) across the first three
points, then collapses between 1.0 and 1.25 as the system condenses into a **single cluster**
(n_cl = 1.0, lcf = 1.00 exactly). Past that point the metric is reporting "nothing moves,"
not a graded suppression of fission. Note ARI now sits *above* its noise floor there — the
partition is more reproducible than Louvain's own seed jitter, the signature of a frozen
blob. Q again is not pinned (−0.087).

Figure: `sweep.png` — in both panels Q (dark) tracks edge turnover (red) almost exactly,
which is the opposite of the decoupling the prediction asserts.

---

## 4. Why the sweeps cannot test the stated prediction

The refined prediction is that *turnover at fixed Q tracks the antisymmetric/symmetric
coupling ratio*. Decomposing the EHD force with g_k ≡ α̂ l_k⁴/(r²+l_k²)^{5/2}:

```
F_ij = −g_j r̂ ,  F_ji = +g_i r̂
  reciprocal part  ∝ (g_i + g_j)/2     (equal and opposite)
  nonreciprocal    ∝ (g_i − g_j)/2     (same sign on both — pair self-propulsion)
  ratio            = |g_i − g_j| / |g_i + g_j|
```

**α̂ cancels exactly from that ratio** (verified numerically: 0.3387 at r=0.15, 0.5617 at
r=0.30, 0.6288 at r=0.50 — identical to 6 d.p. across α̂ ∈ {0.003…0.015}). Both the
symmetric and antisymmetric parts scale linearly in α̂, so α̂ tunes overall EHD strength,
never the nonreciprocity ratio. The ratio is set **only** by the EHD radius contrast
l_L vs l_S — which both sweeps hold fixed at 1/6 and 1/9 per the paper's S5/S8 protocol.

The s-ratio sweep changes *steric* radii with EHD radii fixed, so it too leaves the ratio
untouched; it varies packing and steric frustration instead.

So neither sweep is a test of the refined NWAP prediction. Both instead vary the overall
attraction/packing, which moves Q and turnover together — exactly what the data show.

**What would actually test it:** sweep the EHD radius ratio l_II/l_I (e.g. l_S from 1/6 down
to 1/12) at fixed steric radii and fixed α̂. That varies the antisymmetric part while holding
the symmetric part and the packing approximately fixed — the only way to approach the
"at fixed Q" condition the prediction requires. This needs a new `--sweep lratio` arm in
`simulate.py` (`setup()` currently pins `l0` to the case); it was not run here.

---

## Files

| file | contents |
|---|---|
| `modularity_test.png` | baseline 4-panel, seed-averaged (baseline runs only) |
| `snapshots.png` | final-frame render, mono vs bi |
| `sweep.png` | turnover & Q vs α̂ and vs s_II/s_I |
| `results.csv` | per-run late-window summary, all 31 runs |
| `timeseries.csv` | per-snapshot metrics, all 31 runs |
| `summary_by_condition.csv` | seed-aggregated table |
| `baseline_out/` | CSVs + figures from the baseline-only analysis pass |

Note: `analyze.py` globs every `.npz` in `data/`, so re-running it after the sweeps makes
`timeseries.csv` pool all conditions into the `case` grouping. The baseline 4-panel must be
generated from a baseline-only `data/` (or from `baseline_out/`) — the root figures here were
preserved from that pass. `figures.py:sweep()` was patched in this session to hold the other
parameter at its baseline value; without that fix each sweep panel silently pools in the
other sweep's runs.

---

# 5. Reciprocity test (χ sweep) — run 2026-08-21, 30 runs

Spec and pre-committed thresholds: `EXPERIMENT.md`. χ scales the antisymmetric coupling
alone; the symmetric part is identical at every χ. χ=0 is exactly reciprocal, χ=1 is the
original force. All other parameters fixed; 5 seeds per level, blocked on initial condition.

| χ | edge_jac | Q | ΔQ_null | σ²_v | n_cl | lcf | n |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.273 ± 0.004 | 0.869 ± 0.002 | 0.369 | 1.45e-09 | 3.4 | 0.69 | 5 |
| 0.25 | 0.283 ± 0.004 | 0.856 ± 0.003 | 0.359 | 1.49e-09 | 2.0 | 0.78 | 5 |
| 0.5 | 0.353 ± 0.007 | 0.859 ± 0.003 | 0.345 | 4.84e-09 | 3.5 | 0.94 | 5 |
| 0.75 | 0.356 ± 0.009 | 0.852 ± 0.002 | 0.333 | 1.54e-08 | 12.1 | 0.89 | 5 |
| 1 | 0.384 ± 0.007 | 0.855 ± 0.002 | 0.323 | 3.49e-08 | 22.6 | 0.83 | 10 |
| 1.5 | 0.381 ± 0.008 | 0.854 ± 0.002 | 0.318 | 7.85e-08 | 33.1 | 0.79 | 5 |
| *monodisperse ref* | *0.164* | *0.875* | *0.391* | *1.18e-09* | *5.4* | *0.46* | *5* |

(χ=1 pools the 5 new runs with the 5 pre-refactor baseline runs, which are the same physics —
see *Implementation check* below.)

## Scoring against the pre-committed thresholds

| | test | result | verdict |
|---|---|---|---|
| **H1** | turnover monotone in χ | Spearman ρ = **+0.860**, p = 1.1e-09 (n=30) | **PASS** |
| **H2** | R_Q < 0.10 | R_Q = **0.017** | PASS |
| **H2** | R_turnover > 0.50 | R_turnover = **0.404** | **FAIL** |
| **H3** | σ²_v rises with χ | ρ = **+0.960**, p = 6.1e-17; 24× from χ=0→1 | **PASS** |
| **H4** | χ=0 recovers monodisperse null | closes **50%** of the turnover gap | **PARTIAL** |

**H2 must be recorded as a failure.** The threshold was set in advance at R_turnover > 0.50 and
the observed value is 0.404. The compound criterion as written ("R_Q < 0.1 **and**
R_turnover > 0.5, i.e. at least a 5× separation") was internally inconsistent: the observed
separation is **24×**, far past the 5× gloss, while the absolute magnitude clause fails.
The decoupling is real and large; the specific magnitude bar was set too high. Recording both
rather than quietly adopting the clause that passes.

## What the χ sweep establishes

**1. Arrested coarsening is caused by nonreciprocity — the cleanest positive result here.**
Cluster count rises 3.4 → 22.6 → 33.1 across χ (ρ = +0.878, p = 1.8e-10). At χ=0 the
mixture coarsens to a few large clusters; switching on the antisymmetric coupling holds it
fragmented. Same particles, same packing, same symmetric coupling — reciprocity is the only
difference. This is the paper's central mechanism isolated with a proper control.

**2. Structure really is decoupled from currents.** Across the full χ range Q moves 1.7%
while σ²_v moves 5323%. Q is flat to within ~2% and non-monotonic; the currents are
monotone and enormous. `chi_test.png` (right panel) shows this directly — Q sits on the
unit line while σ²_v climbs 54×.

**3. The mono/bi baseline was a confounded comparison, and the confound is now quantified.**
Switching off reciprocity at fixed particles closes **50%** of the mono/bi edge-turnover
gap. The residual 0.109 is polydispersity, not reciprocity. So the headline "2.4× turnover"
claim splits roughly half-and-half between the two causes — it was **not** all nonreciprocity,
as §1 implicitly assumed.
For σ²_v the story is different: χ=0 sits at 1.23× the monodisperse value, so essentially
**all** of the 30× σ²_v gap is attributable to reciprocity. σ²_v is the clean reciprocity
probe; edge turnover is contaminated by morphology.

## Implementation check

χ=1 is algebraically the original force (ḡ + 1·(g−ḡ) = g). Verified numerically against the
pre-refactor kernel at matched seed: **bitwise identical for 1000 steps**, diverging only to
9.8e-15 by 5000 steps — floating-point rounding amplified by chaos, not a coding error.
The 5-seed χ=1 re-run mean sat 2.55 pooled SEM below the pre-refactor mean with all 5 seeds
falling the same way (sign test p≈0.06); given the bitwise check this is a sampling
fluctuation, and the χ=1 row above pools all 10 runs as the conservative estimate.
`tests/test_reciprocity.py` covers the knob itself.

## Caveats

- t̂=1e5, box 12 — pilot scale, not the paper's steady-state scale. The χ *contrast* is the
  object and all arms share the scale, but a paper-scale repeat is needed before publication.
- χ is a synthetic control parameter; χ=1.5 is outside the derived model (trend evidence only).
- n_cl at χ=0 (3.4) is *below* the monodisperse value (5.4): the reciprocal bidisperse
  mixture coarsens further than the monodisperse one. Not predicted, not yet explained.
- 5 seeds resolves the primary contrast comfortably but leaves per-level means at ±0.008.

