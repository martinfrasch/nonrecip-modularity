# Nonreciprocal colloids × Network-Weighted Action: contact-network modularity test

Reimplementation of the agent-based model from **Hara et al., "Arrested coarsening in
active colloidal suspensions driven by nonreciprocal electrohydrodynamic interactions,"
PRL 137, 068302 (2026)** — DOI 10.1103/96ky-d1p9, arXiv:2509.23164 (End Matter Sec. II,
Eqs. 7–9) — plus a network-analysis pipeline testing the NWAP directed-graph prediction.

## Finding from the pilot run (container, 2026-08-09)

Static Newman modularity Q does **not** discriminate reciprocal from nonreciprocal
dynamics (Q ≈ 0.88 in both; excess over degree-preserving null is *larger* for the
frozen monodisperse crystal). The discriminating observables are **currents**:

| metric (late window)             | monodisperse | bidisperse | ratio |
|----------------------------------|-------------:|-----------:|------:|
| Q (Louvain)                      | 0.883        | 0.878      | ~1    |
| ΔQ vs degree-preserving null     | 0.39         | 0.31       | 0.8   |
| edge Jaccard turnover / 9 s      | 0.18         | 0.47       | 2.6   |
| consecutive-partition ARI        | 0.73 (≈noise floor 0.80) | 0.54 (floor 0.76) | — |
| n_clusters (arrested coarsening) | ~5, falling  | ~21, flat  | 4     |
| σ²_v (paper's order parameter)   | 1.4e-9       | 3.7e-8     | 26    |

Interpretation: the antisymmetric (nonreciprocal) part of the coupling is solenoidal —
it cannot alter a static structural observable, so its entire signature appears in
probability currents: **frozen Q, circulating partition**. Refined NWAP prediction:
turnover rate at fixed Q tracks the antisymmetric/symmetric coupling ratio
(→ `--sweep alpha` and `--sweep sizeratio`).

## Setup

```bash
python -m venv venv && source venv/bin/activate     # or conda
pip install -r requirements.txt
```

## Reproduce the pilot, properly (multi-seed)

```bash
python simulate.py --sweep baseline --T 1e5 --workers 8   # 5 seeds × 2 cases
python analyze.py                                          # -> results.csv, timeseries.csv, summary_by_condition.csv
python figures.py                                          # -> modularity_test.png, snapshots.png
```

Runtime guide (single core ≈ 700 s per 10⁶ steps at N=1000; scales ~linearly in N,
jobs parallelize across seeds): baseline sweep ≈ 30–60 min on 8 cores at default size.

## Paper-scale / steady-state runs

The paper's cluster statistics use t̂ > 4×10⁵ after annealing (S5/S8):

```bash
python simulate.py --sweep baseline --N 4000 --box 24 --T 4e5 --nsnap 200 --workers 8
```

## Sweeps (the actual NWAP test)

```bash
python simulate.py --sweep alpha --workers 8       # α̂ ∈ {0.003, 0.005, 0.010, 0.015}
python simulate.py --sweep sizeratio --workers 8   # s_II/s_I: 0.667 (head-large) → 1.5 (tail-large)
python analyze.py && python figures.py --sweep     # -> sweep.png
```

Predicted outcome: edge/partition turnover increases with α̂ and decreases toward the
tail-large regime (cf. paper Fig. 4b, S5: tail-large suppresses fission), while Q stays
pinned near its energy-set value.

## Model/analysis notes & caveats

- Euler–Maruyama, default dt̂ = 0.05 (pilot used 0.12; stability bound dt̂ < 2/(k·max mobility) ≈ 0.2).
- Nondimensional units: length λ = 9 µm, time τ = 0.01 s. σ̂ = 2.6e-3, α̂ = 0.005 baseline.
- Monodisperse null uses N/2 large particles ⇒ matched area packing (~30%).
- Contact edge: r < 1.1 (s_i + s_j), periodic minimum image.
- Community detection: Louvain (seeded). Always report the **same-graph two-seed ARI
  noise floor** next to consecutive-snapshot ARI; the edge-Jaccard turnover is
  algorithm-free and is the headline current metric.
- Known model limitation (theirs, inherited): pairwise forces only — no many-body
  hydrodynamics, which shifts absolute cluster sizes (paper SI Sec. S6).
- Pilot ran one seed per condition, 12×12 box, t̂ = 5×10⁴; treat pilot numbers as
  provisional until the baseline sweep confirms them.

## Files

- `simulate.py` — numba model, CLI, multi-seed/multi-process, sweeps, resumable output dir
- `analyze.py`  — contact graphs, Q, rewired null, ARI + noise floor, edge turnover, σ²_v
- `figures.py`  — seed-averaged 4-panel figure, snapshot render, sweep figure
- `pilot/`      — raw pilot outputs from the container session (figures + results.json)
