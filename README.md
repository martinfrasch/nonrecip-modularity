# Isolating the antisymmetric sector of a nonreciprocal colloidal model

Code, data summaries and manuscript for

> M. G. Frasch, *Isolating the antisymmetric sector of a nonreciprocal colloidal model: kinetic
> unjamming, transient arrested coarsening, and irreversibility without coarse-grained
> circulation* (2026). Source: `PAPER.md`; typeset: `paper.pdf`; arXiv package:
> `arxiv-submission.tar.gz`.

The model is a reimplementation of the agent-based colloid model of Hara et al., *Phys. Rev.
Lett.* **137**, 068302 (2026), arXiv:2509.23164 (End Matter Sec. II, Eqs. 7–9): overdamped
size-asymmetric colloids driven by electrohydrodynamic flows, with a nonreciprocal pair force. We
add a **reciprocity mixing parameter χ** that scales the antisymmetric part of the pair coupling
while leaving the symmetric part unchanged (χ = 0 restores Newton's third law exactly, χ = 1
recovers the published law), and a set of network, kinetic and thermodynamic measurements built
on it.

## What the paper reports

- Scaling the antisymmetric sector alone changes equal-time structure by an order of magnitude;
  it unjams the reciprocal gel at weak drive and sustains a fragment population at strong drive,
  so cluster count is non-monotonic in χ.
- The fragmented state is a long-lived transient: the largest cluster grows as N^1.00 ± 0.07 and
  the condensation time as N^1.0 over a fourfold range of N.
- The selected scale, S* ≈ 8–10 particles, is a critical-size crossover: the committor of the
  fragment-size dynamics, converged in the horizon, is one half at S ≈ 8, and the
  monomer-exchange drift changes sign at S ≈ 10.
- Excess dissipation over the reciprocal control, from a configurational estimator that needs no
  trajectory, is unresolved at χ ≤ 0.25 and quadratic from χ = 0.75; the symmetric forces cancel
  87–91 % of the unopposed antisymmetric contribution above that crossover.
- The fragment population is a kinetic balance whose channel pairs (shedding/reabsorption,
  fission/fusion, nucleation/dissolution) are each approximately balanced.
- On single-particle coordinates a single effective temperature equal to the bath temperature is
  shared by both species over a short-lag window that shrinks with drive.
- No circulation is detectable in any coarse observable examined, including the cluster-size
  coordinate and an all-pairs lag test; the antisymmetric pair law is predominantly gradient on
  the contact graph. Newman modularity is degenerate on these networks.

## Setup

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

Runtime ≈ 700 s per 10⁶ steps at N = 1000 on one core, scaling linearly in N; jobs parallelise
across seeds. Raw trajectories (`data*/`, ~4 GB) are not version-controlled.

## Pipeline

```bash
# simulate -> data/*.npz (resumable; existing files are skipped)
python simulate.py --sweep chi --workers 8                          # pilot scale, N = 1000
python simulate.py --sweep chipaper --workers 8                     # paper scale, N = 4000
python continue_run.py --pattern "data_paper/*_x1.5_*.npz" --T 4e5  # continuations
python continue_run.py --pattern "..." --T 2.5e4 --nsnap 500 --outdir data_dense  # dense kinetics

# network observables -> results.csv, timeseries.csv, summary_by_condition.csv
python analyze.py --burn 0.5

# kinetics, fates, dissipation, response
python static_epr.py --pattern "data_paper/*.npz" --out epr_static.csv \
    --calib-dt 0.005 --calib-n 2 --calib-T 500 --calib-nsamp 100      # configurational dissipation
python cluster_kinetics.py --pattern "data_dense/*.npz"               # event classification, drift
python fragment_budget.py --pattern "data_dense/*seed1*.npz"          # source-sink budget
python committor.py --pattern "data_paper/*cont2_dense500.npz" --T 3000
python hs_fdt.py --pattern "data_paper/*_x1_N4000*seed[123].npz" && python hs_analysis.py
python coarsening_law.py --pattern "data_paper/*.npz" "data_box/*.npz"
python hodge_contact.py --pattern "data_paper/*_N4000_*seed[123].npz"
python lag_asymmetry.py --pattern "data_dense/*.npz"

# figures and manuscript
python make_figures.py                     # fig1..fig6
python build_tex.py && latexmk -pdf paper.tex
```

## Files

Simulation: `simulate.py` (model, χ, sweeps), `simulate_par.py` (thread-parallel kernel),
`simulate_rect.py` (rectangular box), `continue_run.py`. Analysis: `analyze.py`, `static_epr.py`,
`epr_probe.py`, `cluster_kinetics.py`, `fragment_budget.py`, `committor.py`, `hs_fdt.py`,
`hs_analysis.py`, `fdt.py`, `fdt_driftfree.py`, `drift_control.py`, `onsager*.py`,
`coarsening_law.py`, `hodge_contact.py`, `lag_asymmetry.py`, `crossover.py`, `arrest_test.py`,
`irreversibility.py`, `benchmark_irrev.py`. Comparison systems: `cilia_analysis.py`,
`single_cluster_circulation.py`, `vicsek.py`, `hodge_dominance.py`, `forkosh_asymmetry.py`.
Figures: `make_figures.py`, `figures.py`. Manuscript: `PAPER.md`, `build_tex.py`, `paper.tex`.
Validation: `tests/test_reciprocity.py`.

Per-run summaries are the `*.csv` files at the top level; `pilot/` and `baseline_out/` hold the
single-seed pilot and the baseline sweep. `EXPERIMENT.md` is the pre-registered protocol with
thresholds committed before any χ run; `MANUSCRIPT_AUDIT.md` records the internal review rounds
and every correction; the remaining `*.md` files are working notes and per-analysis result
records referenced from the manuscript's data-availability section.

Model conventions that constrain edits are in `CLAUDE.md`.
