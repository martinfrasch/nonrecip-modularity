# Isolating the antisymmetric sector of a nonreciprocal colloidal model

Code, data summaries and manuscript for

> M. G. Frasch, *Isolating the antisymmetric sector of a nonreciprocal colloidal model: kinetic
> unjamming, transient arrested coarsening, and irreversibility without coarse-grained
> circulation* (2026). Source `paper/PAPER.md`, typeset `paper/paper.pdf`, arXiv package
> `paper/arxiv-submission.tar.gz`.

## Why this study exists

[minAction.net](https://minaction.net) is a programme asking whether living organisation —
physiological regulation, the modular structure of metabolic networks, learned models of physical
systems — can be read as minimising a **Network-Weighted Action**,

    S_NW = ∫ (E − I + A·C) dt ,

an energy-like cost E traded against information I under the cost C of forming and maintaining
connections, weighted by connectivity A (Frasch 2026, *J. Physiol.*, DOI 10.1113/JP290762). Its
published tests so far are meta-analytic or computational. A framework of that kind needs a
physical system in which one knob varies the quantity it singles out and a prediction can fail
cleanly.

Nonreciprocal colloids are that system. In the agent-based model of Hara et al. (*Phys. Rev. Lett.*
**137**, 068302, 2026; arXiv:2509.23164) size-asymmetric colloids driven by electrohydrodynamic
flows exert unequal forces on each other. We add a **reciprocity mixing parameter χ** that scales
the antisymmetric part of the pair coupling while leaving the symmetric part bit-for-bit unchanged:
χ = 0 restores Newton's third law exactly, χ = 1 recovers the published force law. Everything else
— particles, packing, noise stream, symmetric attraction — is held fixed.

## The hypotheses, and what happened to them

The hypothesis under test is a *directed-graph extension* of NWAP, formulated by the author while
designing this study and not published elsewhere: split a nonreciprocal coupling into a symmetric
sector that selects structure by energy minimisation and an antisymmetric sector that adds
circulation on top without disturbing it. Its consequences, and four pre-registered predictions
(`paper/EXPERIMENT.md`, thresholds committed before any χ run), were tested against the outcome.

| hypothesis | prediction | outcome |
|---|---|---|
| **Solenoidality** — the antisymmetric sector is structure-preserving | no equal-time structural observable moves with χ | **Refuted.** Cluster count changes 12× with the symmetric coupling fixed; the antisymmetric drift does not preserve the reciprocal stationary density (a pair calculation, not only a simulation). |
| **Circulating partition** — the sector's signature is in probability currents | signed-area rates in coarse observable planes, cluster-size cycles, lag asymmetries | **Refuted, in every projection examined.** The dynamics are irreversible (resolved excess dissipation) but nothing coarse circulates; the antisymmetric pair law is itself predominantly gradient on the contact graph. |
| **Modularity excess** (Frasch 2026, arXiv:2605.05254) — constrained organisation shows as sustained Newman-modularity excess over a null | excess grows with nonreciprocity | **Refuted with inverted sign.** Modularity is degenerate on these contact networks and its excess is largest in the reciprocal reference. |
| **Structure set by the antisymmetric/symmetric ratio** | cluster statistics controlled by χ | **Confirmed in substance, non-monotonically.** Weak nonreciprocity unjams the reciprocal gel (fivefold *fewer* clusters at χ = 0.25), strong nonreciprocity fragments it (two orders of magnitude more); no proposed functional predicts the dip, a balance of kinetic rates describes it. |
| **Near-equilibrium (tier-II) description** — quadratic dissipation, Onsager reciprocity, a single effective temperature | all three signatures | **Undecided as a tier; decided as a window.** Quadratic dissipation is near-tautological at fixed structure; an apparent Onsager antisymmetry does not survive a structurally dissimilar control; a single effective temperature equal to the bath's holds on single-particle coordinates over a short-lag window that shrinks with drive. |
| H1 (pre-registered): edge turnover rises monotonically in χ | Spearman ρ > 0 | **Failed at paper scale** (non-monotonic; passed only in the unequilibrated pilot). |
| H2 (pre-registered): "frozen Q, circulating partition" | R_Q < 0.1 and R_turnover > 0.5 | **Failed informatively** — Q stayed flat because it is degenerate, not because structure is preserved. |
| H3 (pre-registered): σ²_v rises with χ | directional | **Passed.** |
| H4 (pre-registered): χ = 0 turnover falls to the monodisperse value | if not, the mono/bi gap was polydispersity | **Partial** — about half the gap is polydispersity; edge turnover retired as a probe. |

What survived is the instrument, not the mechanism: constructing the symmetric/antisymmetric split
is what made every measurement possible, and the sector it isolates has real, measurable
consequences — unjamming, fragment turnover, a critical-size crossover at 8–10 particles,
extensive condensation whose time scales as N, and dissipation that switches on across the same
structural crossover at which the symmetric forces stop cancelling it. None of it validates the
NWAP functional. The structure the sector selects is a kinetic balance — shedding against
reabsorption, fission against fusion, nucleation against dissolution, each approximately balanced —
and that balance is what a functional for structure selection would have to reproduce.

## Next steps

- **Structure selection.** Find, or rule out, a functional whose extremum reproduces the measured
  rate balance. Nonequilibrium quasipotentials and large-deviation constructions are the
  candidates; nothing tested here covers it.
- **Frequency-resolved response.** Compute the correlation and dissipative response spectra and
  the Harada–Sasa sum rule, to turn the short-lag agreement window into a spectral statement and
  tie the fluctuation–dissipation violation to the measured dissipation. Green–Kubo at χ = 0 for
  the low-drive curvature.
- **Hydrodynamics.** The model is pairwise. Many-body hydrodynamics are the ingredient that could
  stabilise finite clusters; a far-field mobility coupling is the minimal test of whether the
  source experiment's arrest is hydrodynamic.
- **Thermodynamic limit.** A fifth system size (N ≈ 40 000) to twice the current duration to
  extend the t_x ∝ N measurement.
- **The physical knob.** Sweep the electrohydrodynamic radius contrast, the experimentally
  accessible analogue of χ, with per-condition compensation of the symmetric coupling.
- **Dominance data.** Raw interaction tracking for the social-dominance analogue, where the
  published matrices are underpowered (`notes/REQUEST_RAW_TRACKING.md`).

## Repository layout

```
sim/        model and integrators: simulate.py, simulate_par.py, simulate_rect.py, continue_run.py
analysis/   measurements: analyze.py, static_epr.py, cluster_kinetics.py, fragment_budget.py,
            committor.py, hs_fdt.py, hs_analysis.py, coarsening_law.py, hodge_contact.py,
            lag_asymmetry.py, fdt*.py, onsager*.py, epr_probe.py, crossover.py, ...
external/   comparison systems and their open data: cilia_analysis.py, vicsek.py,
            hodge_dominance.py, forkosh_asymmetry.py, data_dominance/, data_forkosh/
figures/    make_figures.py and the manuscript figures
paper/      PAPER.md (source of record), build_tex.py, paper.tex, paper.pdf, arXiv package,
            EXPERIMENT.md (pre-registration), MANUSCRIPT_AUDIT.md (review rounds and corrections)
results/    per-run summaries (CSV), pilot/ and baseline_out/
notes/      working notes and per-analysis result records
tests/      test_reciprocity.py (the χ construction)
```

## Setup and pipeline

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python -m tests.test_reciprocity

# simulate -> data*/ (gitignored; resumable, existing files are skipped)
python -m sim.simulate --sweep chi --workers 8                        # pilot scale, N = 1000
python -m sim.simulate --sweep chipaper --workers 8                   # paper scale, N = 4000
python -m sim.continue_run --pattern "data_paper/*_x1.5_*.npz" --T 4e5
python -m sim.continue_run --pattern "data_paper/*_x1_*seed[123].npz" --T 2.5e4 --nsnap 500 --outdir data_dense

# measurements -> results/
python -m analysis.analyze --burn 0.5
python -m analysis.static_epr --pattern "data_paper/*.npz" --calib-dt 0.005 --calib-n 2 --calib-T 500 --calib-nsamp 100
python -m analysis.cluster_kinetics --pattern "data_dense/*.npz"
python -m analysis.fragment_budget --pattern "data_dense/*seed1*.npz"
python -m analysis.committor --pattern "data_paper/*cont2_dense500.npz" --T 3000
python -m analysis.hs_fdt --pattern "data_paper/*_x1_N4000*seed[123].npz" && python -m analysis.hs_analysis
python -m analysis.coarsening_law --pattern "data_paper/*.npz" "data_box/*.npz"
python -m analysis.hodge_contact --pattern "data_paper/*_N4000_*seed[123].npz"
python -m analysis.lag_asymmetry --pattern "data_dense/*.npz"

# figures and manuscript
python -m figures.make_figures
python paper/build_tex.py && (cd paper && latexmk -pdf paper.tex)
```

Runtime ≈ 700 s per 10⁶ steps at N = 1000 on one core, linear in N; jobs parallelise across
seeds. Raw trajectories (~4 GB) are not version-controlled. The cilia dataset is open on Dryad
(10.5061/dryad.0gb5mkm2j, CC0); the dominance matrices and the Forkosh behavioural table ship
in `external/`.
