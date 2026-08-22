# Results — baseline + sweeps, 31 runs, M2 Max (2026-08-21)

Environment: Python 3.12 (arm64), numba 0.67.0, numpy 2.5.2, networkx 3.6.1, pandas 3.0.5.
All runs N=1000 (mono: 500), box 12, t̂=1e5, dt̂=0.05, 100 snapshots, `--burn 0.5`.
31 `.npz` total: baseline 5 seeds × 2 cases, α̂ sweep 4×3, s-ratio sweep 5×3.

## Verdict in one line

The baseline "frozen Q, circulating partition" claim **reproduces cleanly**. The refined
sweep prediction **does not hold as stated** — and the reason is that neither sweep
actually varies the quantity the prediction is about.

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
