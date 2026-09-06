# Nonreciprocal colloids × Network-Weighted Action: contact-network modularity test

Reimplementation of the agent-based model from **Hara et al., "Arrested coarsening in
active colloidal suspensions driven by nonreciprocal electrohydrodynamic interactions,"
PRL 137, 068302 (2026)** — DOI 10.1103/96ky-d1p9, arXiv:2509.23164 (End Matter Sec. II,
Eqs. 7–9) — plus a network-analysis pipeline built to test the NWAP directed-graph
prediction. **Note:** the pipeline as written builds *undirected* contact graphs
(`nx.Graph`) throughout, so the directed content of that prediction remains untested.

## Results (113 independent trajectories + 25 continuations, N=1000 to N=22,000; see `PAPER.md`, `FINAL_RESULTS.md`)

**Nonreciprocity causes arrested coarsening, and that is a change in *static* structure.**
Scaling the antisymmetric coupling alone (χ, see below) with the reciprocal part held
bit-identical takes the system from 11 clusters to 137 at paper scale — 12.4× —
and raises activity σ²_v by 133×. This reproduces the paper's central
mechanism with a control the paper does not have: a same-particle, same-packing,
same-symmetric-coupling reciprocal reference.

| paper scale, 3 seeds | χ=0 (reciprocal) | χ=1 (original force) | χ=1.5 |
|---|---:|---:|---:|
| n_clusters | 11.1 | 95.2 | 137.3 |
| σ²_v | 4.19e-10 | 2.76e-08 | 5.57e-08 |
| Newman Q | 0.930 | 0.906 | 0.901 |
| edge turnover | 0.260 | 0.378 | 0.349 |

### The original interpretation did not survive

The pilot's reading — *"the antisymmetric part is solenoidal, it cannot alter a static
structural observable, so its entire signature appears in probability currents: frozen Q,
circulating partition"* — fails on both halves:

- **Q is frozen because Q is degenerate here, not because structure is invariant.** Across
  all conditions run, cluster count spans 40× while Q spans 14%, and the two are
  uncorrelated (ρ=+0.31, p=0.28). A 3-cluster and a 33-cluster state are indistinguishable
  in Q. Louvain's resolution limit merges anything smaller than ~√(2E) nodes, and this
  holds at 4× the system size.
- **Static structure does change**, by 12.4× in `n_cl` — a single-snapshot observable.
  It was already visible in the pilot table below as `n_clusters ~5 vs ~21`.
- **No circulation has been detected.** The signed area rate in observable planes is exactly
  zero under detailed balance; measured at both scales, every plane at every χ is
  indistinguishable from the equilibrium null (all p ≥ 0.38).
- **The system *is* genuinely irreversible.** Entropy production — which unlike edge turnover
  and σ²_v is exactly zero under detailed balance — is positive for χ ≳ 0.75 and scales
  quadratically in χ at high drive. So a current exists; it simply does not appear as
  circulation in any network observable measured here. See `PAPER.md` §3.10 for what this does
  and does not establish about near-equilibrium response — in short, less than we first thought.

### Why the original sweeps could not test the refined prediction

α̂ **cancels exactly** from the antisymmetric/symmetric coupling ratio: both parts scale
linearly in it, so α̂ tunes overall EHD strength, not reciprocity. The s_II/s_I sweep moves
steric radii with the EHD radii pinned, so it varies packing instead. Both sweeps move
structure and currents together, and neither varies the quantity the prediction is about.

The knob that does is **χ**, a reciprocity mixing parameter that scales the antisymmetric
part while leaving the symmetric part identically unchanged (χ=0 exactly reciprocal, χ=1 the
original Hara et al. force). See `EXPERIMENT.md` for its spec and `tests/test_reciprocity.py`
for its validation.

### Metric health

`n_cl` and σ²_v are the robust observables. **Edge Jaccard turnover is not**: it correlates
with morphology (ρ=+0.77 pilot, +0.62 paper scale, against n_cl), 52% of its mono-vs-bi gap is
polydispersity rather than reciprocity, ~70% of its value survives in the provably current-free
χ=0 system, and it turns
**non-monotonic** in χ at paper scale. Newman Q should not be used as a discriminator on
these graphs at all. `AUDIT.md` documents this in full.

<details>
<summary>Original pilot table (2026-08-09, single seed) — retained for provenance</summary>

| metric (late window)             | monodisperse | bidisperse | ratio |
|----------------------------------|-------------:|-----------:|------:|
| Q (Louvain)                      | 0.883        | 0.878      | ~1    |
| ΔQ vs degree-preserving null     | 0.39         | 0.31       | 0.8   |
| edge Jaccard turnover / 9 s      | 0.18         | 0.47       | 2.6   |
| consecutive-partition ARI        | 0.73 (≈noise floor 0.80) | 0.54 (floor 0.76) | — |
| n_clusters (arrested coarsening) | ~5, falling  | ~21, flat  | 4     |
| σ²_v (paper's order parameter)   | 1.4e-9       | 3.7e-8     | 26    |

Note this comparison confounds reciprocity with polydispersity: monodisperse and bidisperse
differ in both. The χ=0 arm is the control that separates them.

</details>

## Manuscript

`PAPER.md` is the source of record. `build_tex.py` converts it to REVTeX 4.2 (APS/PRE) and the
PDF is built with `latexmk`:

```bash
python build_tex.py && latexmk -pdf paper.tex     # -> paper.tex, paper.pdf
```

The conversion is mechanical: the markdown carries physics notation in Unicode rather than LaTeX,
so `build_tex.py` protects code spans and citations, maps Unicode to math, escapes what remains,
and emits the 20 tables as captioned floats and the 6 figures at the sections that discuss them.
Equation blocks are hand-mapped (`EQ` in the script) and table captions live in `CAPS`.
`arxiv-submission.tar.gz` is the arXiv package: single self-contained `.tex` with an inline
bibliography plus the six PNGs, no `.bbl` needed.

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

## Sweeps (superseded — see `--sweep chi` below)

```bash
python simulate.py --sweep alpha --workers 8       # α̂ ∈ {0.003, 0.005, 0.010, 0.015}
python simulate.py --sweep sizeratio --workers 8   # s_II/s_I: 0.667 (head-large) → 1.5 (tail-large)
python analyze.py && python figures.py --sweep     # -> sweep.png
```

These two sweeps were run and are **retained for the record only**: α̂ cancels from the
antisymmetric/symmetric ratio and s_II/s_I leaves the EHD radii fixed, so neither varies
reciprocity (RESULTS.md §4). Observed: turnover *falls* with α̂ (0.601 → 0.282, opposite to
the original prediction) and Q is not pinned in either sweep. Use `--sweep chi` instead.

```bash
python simulate.py --sweep chi --workers 8      # the actual reciprocity test
python analyze.py && python figures.py --chi    # -> chi_test.png
python epr_probe.py                             # entropy production (the current probe)
```

## Model/analysis notes & caveats

- Euler–Maruyama, default dt̂ = 0.05 (pilot used 0.12; stability bound dt̂ < 2/(k·max mobility) ≈ 0.2).
- Nondimensional units: length λ = 9 µm, time τ = 0.01 s. σ̂ = 2.6e-3, α̂ = 0.005 baseline.
- Monodisperse null uses N/2 large particles ⇒ matched area packing (~30%).
- Contact edge: r < 1.1 (s_i + s_j), periodic minimum image.
- Community detection: Louvain (seeded). Newman Q is **degenerate on these graphs** and
  should not be used to discriminate conditions; its resolution limit merges anything below
  ~√(2E) nodes. Report the same-graph two-seed ARI noise floor next to any ARI, and note the
  floor is itself condition-dependent (0.35–0.87 observed).
- Edge Jaccard turnover is algorithm-free but **not** a current: ~70% of it survives at χ=0,
  where detailed balance makes the probability current exactly zero. It is also
  morphology-correlated and non-monotonic in χ at paper scale. Entropy production
  (`epr_probe.py`) is the current probe; σ²_v is the robust activity probe.
- Known model limitation (theirs, inherited): pairwise forces only — no many-body
  hydrodynamics, which shifts absolute cluster sizes (paper SI Sec. S6).
- The pilot (one seed, 12×12 box, t̂ = 5×10⁴) was confirmed at 5 seeds but its *interpretation*
  was not — see Results above.
- t̂ = 1×10⁵ at box 12 is **not** steady state: σ²_v contrast grows from 24× there to 133× at
  t̂ = 4×10⁵, and edge turnover changes from monotonic in χ to non-monotonic. Use paper scale
  for any quantitative claim.
- σ²_v as implemented is the variance of *speed* over one snapshot interval (20,000 steps),
  not an instantaneous velocity variance; it is not directly comparable to Hara et al. Fig. 2d.

## Files

- `simulate.py` — numba model, CLI, multi-seed/multi-process, sweeps, resumable output dir
- `analyze.py`  — contact graphs, Q, rewired null, ARI + noise floor, edge turnover, σ²_v
- `figures.py`  — seed-averaged 4-panel figure, snapshot render, sweep figure
- `epr_probe.py` — entropy production at small dt from equilibrated configs (the current probe)
- `pilot/`      — raw pilot outputs from the container session (figures + results.json)
- `FINAL_RESULTS.md` — headline results, both scales
- `AUDIT.md`    — design critique: which observables test which hypothesis, and which do not
- `EXPERIMENT.md` — χ reciprocity-test spec with pre-committed thresholds
- `RESULTS.md`  — pilot-scale record, including the two mis-specified sweeps
