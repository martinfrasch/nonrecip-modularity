# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A reimplementation of the agent-based colloid model from Hara et al., PRL 137, 068302 (2026)
(arXiv:2509.23164, End Matter Sec. II, Eqs. 7–9), plus a contact-network analysis pipeline
built to test a prediction from the Network-Weighted Action Principle (NWAP) directed-graph
extension. It is a research pipeline, not an application: four files, no package, no tests,
no build step.

Not a git repo yet — `bootstrap_github.sh` does `git init` + `gh repo create` under the
**martinfrasch** account.

## Commands

```bash
python -m venv venv && source venv/bin/activate && pip install -r requirements.txt

# simulate -> data/*.npz  (resumable: existing .npz files are skipped, so reruns are cheap)
python simulate.py --case bidisperse --seed 12          # single run
python simulate.py --sweep baseline --T 1e5 --workers 8 # 5 seeds x 2 cases, ~30-60 min on 8 cores
python simulate.py --sweep alpha --workers 8            # alpha in {0.003,0.005,0.010,0.015}
python simulate.py --sweep sizeratio --workers 8        # s_II/s_I from 0.667 (head-large) to 1.5

# analyze -> results.csv, timeseries.csv, summary_by_condition.csv
python analyze.py --burn 0.5

# figures -> modularity_test.png, snapshots.png (or sweep.png)
python figures.py
python figures.py --sweep
```

Runtime ≈ 700 s per 10⁶ steps at N=1000 on one core; scales ~linearly in N, parallelizes
across seeds. To smoke-test a change, use a tiny run (`--N 200 --box 6 --T 1e3 --nsnap 10`)
rather than a full sweep.

## Pipeline shape

`simulate.py` → `data/*.npz` → `analyze.py` → `*.csv` → `figures.py` → `*.png`

The stages are coupled only through those files. Filenames encode the condition
(`{case}_a{alpha}_sr{s_ratio}_N{N}_L{box}_T{T}_seed{n}.npz`) and each `.npz` carries its own
metadata (`case`, `seed`, `alpha`, `s_ratio`, `Lbox`, `dt_snap`, `svec`, `lvec`, `types`).
`analyze.py` reads that metadata into every CSV row, so if you add a swept parameter you must
touch three places: `build_jobs`, the `tag`/`np.savez_compressed` call in `run_one`, and the
`meta` dict + `groupby` in `analyze.py`.

## Physics that constrains edits

- The nonreciprocity lives in one place: in `_run`, the EHD force on `i` from `j` uses
  `l4[j]` while the force on `j` from `i` uses `l4[i]`. Unequal EHD radii ⇒ action-reaction
  is broken. Symmetrizing that (e.g. reusing one `g` for both) silently destroys the entire
  experiment. The collision spring above it *is* reciprocal (`F[i] += ; F[j] -=`) — that's correct.
- `_run` is a single `@njit` kernel with a hand-rolled cell list; `RCUT = 1.0` is baked into
  the cell sizing *and* the `r2 < 1.0` interaction cutoff. Changing the EHD cutoff means
  changing both.
- Euler–Maruyama; stability bound is `dt < 2/(k·max mobility) ≈ 0.2`. Default `dt = 0.05`.
  The pilot used 0.12 — pilot numbers are provisional.
- The `monodisperse` case runs `N//2` particles by design (matched ~30% area packing against
  the bidisperse mix). Don't "fix" that.
- Nondimensional units throughout: λ = 9 µm, τ = 0.01 s, σ̂ = 2.6e-3, α̂ = 0.005 baseline.
- Periodic box everywhere: minimum-image in `_run`, `cKDTree(boxsize=Lbox)` in
  `contact_graph`, `d -= Lbox*round(d/Lbox)` in the displacement/σ²_v calc.

## Analysis conventions

- Contact edge iff `r < 1.1 * (s_i + s_j)`. The kd-tree prefilter uses `2*svec.max()`, so it
  over-collects and then filters exactly — keep that two-step structure if you change radii.
- **Q is not the discriminator.** The pilot found Q ≈ 0.88 for both reciprocal and
  nonreciprocal dynamics; the antisymmetric coupling is solenoidal and cannot move a static
  structural observable. The headline metric is **edge Jaccard turnover** (algorithm-free);
  Louvain ARI is secondary and must always be reported next to the same-graph two-seed ARI
  noise floor (`ari_floor`), which is ~0.76–0.80 — an ARI of 0.73 is near noise, not a signal.
- `nullQ` (degree-preserving double-edge-swap) is computed only every `NULL_EVERY=4`
  snapshots and is `NaN` otherwise; anything aggregating it must `dropna` or use
  `mean(skipna)`.

## Notes

- `pilot/` holds raw outputs from a single-seed container run (12×12 box, t̂=5×10⁴) kept as a
  reference point. Treat those numbers as provisional; don't overwrite them.
- `data/`, `venv/`, `*.csv`, `*.png` are gitignored by `bootstrap_github.sh` (with
  `!pilot/*.png` excepted).
