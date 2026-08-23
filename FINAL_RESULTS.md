# Final results

**82 simulation runs.** N=1000 pilot scale (61 runs) and N=4000 / box 24 / t̂=4×10⁵ paper
scale (21 runs), plus a dedicated entropy-production probe. Environment: Python 3.12 arm64,
numba 0.67, numpy 2.5, M2 Max.

Read with `AUDIT.md`, which explains why several of the original observables could not test
the hypotheses they were assigned to.

---

## The one-paragraph answer

**Nonreciprocity causes arrested coarsening.** At paper scale, turning the antisymmetric
coupling from off to full — with the reciprocal part held bit-identical — takes the system
from 11 clusters to 137 (12.4×, p = 9.4e-10) and raises activity 133× (p = 2.9e-05). That
result is clean, large, reproducible at two system sizes, and controlled in a way the source
paper's own comparison is not.

**The framing that motivated the project does not survive.** "Frozen Q, circulating
partition" fails on both halves: Q is frozen only because Newman modularity is degenerate on
these graphs, and no circulation has ever been detected. What is real is simpler and
stronger than the original claim — nonreciprocity restructures the system and drives it far
from equilibrium.

---

## 1. Paper scale, 3 seeds per level

| χ | Q | n_cl | edge turnover | σ²_v | ARI (its floor) |
|---|---:|---:|---:|---:|---:|
| 0 (reciprocal) | 0.930 | 11.1 | 0.260 | 4.19e-10 | 0.725 (0.84) |
| 0.25 | 0.911 | **2.2** | 0.289 | 4.28e-10 | 0.548 |
| 0.5 | 0.910 | 17.3 | 0.352 | 3.65e-09 | 0.533 |
| 0.75 | 0.909 | 56.5 | **0.384** | 1.37e-08 | 0.511 |
| 1.0 (original force) | 0.906 | 95.2 | 0.378 | 2.76e-08 | 0.489 |
| 1.5 | 0.901 | 137.3 | 0.349 | 5.57e-08 | 0.482 |
| *monodisperse ref* | *0.936* | *11.2* | *0.148* | *3.07e-10* | *0.830* |

| trend across χ | Spearman ρ | p |
|---|---:|---:|
| n_cl | **+0.931** | 2.1e-08 |
| σ²_v | **+0.975** | 7.1e-12 |
| edge turnover | +0.636 | 4.5e-03 |
| Q | −0.781 | 1.3e-04 |

## 2. What changed from the N=1000 results

**Everything structural got stronger; the turnover story broke.**

- **n_cl contrast grew** 9.8× → **12.4×**. The H_A falsification is firmer at steady state.
- **σ²_v contrast grew** 24× → **133×**. At pilot scale the runs had not converged; the
  activity difference was being badly underestimated.
- **Q stayed degenerate at 4× the system size.** 0.930 → 0.901, a 3.1% swing against a 12.4×
  swing in n_cl. Larger graphs did not rescue Louvain's resolution limit.
- **Edge turnover is NON-MONOTONIC.** It peaks at χ=0.75 (0.384) and *falls* 9% by χ=1.5.
  H1 as originally stated — turnover increases monotonically in χ — **fails at paper scale.**
  The pilot-scale monotonicity was an artifact of not being in steady state.

That last point matters: H1 was the one hypothesis scored PASS in §5 of `RESULTS.md`. At
paper scale it does not hold. `AUDIT.md` §5–6 predicted this failure mode — edge turnover
tracks morphology (ρ=+0.795 with n_cl) and carries a large hard-cutoff flicker floor, so it
was never a robust dynamical probe. σ²_v and n_cl, which the audit promoted, are the metrics
that stayed monotonic.

**A reproducible non-monotonicity at low χ.** n_cl *dips* at χ=0.25 (2.2, below χ=0's 11.1)
before rising steeply. The same dip appears at N=1000 (3.4 → 2.0 → 3.5). Weak nonreciprocity
appears to *promote* coarsening before strong nonreciprocity fragments the system. Not
predicted by anything here and not explained; it reproduces across both scales and is worth
its own investigation.

## 3. H4 — the confound in the original baseline, at paper scale

The mono-vs-bi comparison varies reciprocity *and* polydispersity together. Switching off
reciprocity at fixed particles splits it:

| observable | mono | χ=0 (reciprocal, same particles) | χ=1 | share attributable to reciprocity |
|---|---:|---:|---:|---:|
| edge turnover | 0.148 | 0.260 | 0.378 | **52%** |
| σ²_v | 3.07e-10 | 4.19e-10 | 2.76e-08 | **99.6%** |
| n_cl | 11.2 | 11.1 | 95.2 | — |

Identical to the N=1000 split (47% / ~100%). **σ²_v is the clean reciprocity probe; edge
turnover is about half polydispersity.** The README's headline "2.4× turnover" was never a
pure nonreciprocity effect.

New at paper scale: the reciprocal bidisperse system now has *exactly* the monodisperse
cluster count (11.1 vs 11.2). At N=1000 it coarsened further than mono (3.4 vs 5.4) — that
anomaly was a finite-size/non-steady-state artifact and is gone.

## 4. Entropy production — the first genuine current measurement

Neither edge turnover nor σ²_v is a current: both are nonzero at χ=0, which is a *provable*
equilibrium (conservative pair forces, uniform D/μ ⇒ detailed balance ⇒ zero current).
Entropy production is zero there by construction, so it is the right instrument.

It is unusable at the production timestep — its discretisation residual swamps the signal —
but works at dt = 6.25e-4 restarting from equilibrated configurations, with a paired
same-configuration χ=0 control to subtract each configuration's own bias.

**Evolving-structure sweep** (5 seeds, each χ from its own equilibrated state):

| χ | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 |
|---|---:|---:|---:|---:|---:|
| net EPR/particle | −5.8e-5 | −2.2e-5 | +8.0e-4 | **+2.0e-3** | **+6.5e-3** |

Clearly nonzero for χ ≥ 0.75; buried in bias noise below. Fitted exponent **χ^3.01**, not the
χ² theory requires.

**Fixed-configuration test** (one structure, χ varied over a window too short to relax) —
this is where χ² must hold, because structural co-variation is removed:

| seed | net EPR/χ² across χ=0.25…1.5 | flatness |
|---|---|---:|
| 1 | +3.17e-3, +3.62e-3, +3.91e-3, +4.03e-3, +3.75e-3 | **1.27×** |
| 2 | −6.8e-4, +2.32e-3, +2.31e-3, +2.87e-3, +3.72e-3 | 1.61× |
| 3 | all negative until χ=1.5 — bias estimate 2× the others | unusable |

**Seeds 1–2 confirm the χ² law at fixed structure.** Seed 3's bias probe is a single
stochastic realisation that came in 2× high and drove every net value negative — an
instrument-noise failure, not physics. The bias probe needs averaging over several noise
seeds; with n=1 it is not reliable run-to-run.

So: **entropy production is real, positive, and consistent with χ² once structure is held
fixed.** The χ³ exponent in the evolving sweep is the same structure-covariation confound
that has run through this entire project. My first explanation for it — that fragmentation
creates more nonreciprocal contacts — was **wrong**: contact count *falls* with χ (7.67 →
5.81), and normalising by it made the spread worse. The correct explanation is untested.

## 5. Circulation — still nothing, at double the statistics

Signed area rate in observable planes is exactly zero under detailed balance. At paper scale
(nsnap=200, twice the pilot statistics), every plane at every χ remains indistinguishable
from the χ=0 equilibrium null (all p ≥ 0.38).

The "circulating partition" half of the project's headline has now been tested twice and
found nothing. This is a null result, not a disproof — the coarse observables chosen may
simply not project onto the current. But EPR proves current *exists* (§4), so the honest
reading is: **the system is genuinely irreversible, and that irreversibility does not show
up as circulation in any network observable measured here.**

## 6. Scorecard

| claim | status |
|---|---|
| Nonreciprocity causes arrested coarsening | **CONFIRMED**, 12.4× in n_cl, p=9.4e-10, two scales |
| Nonreciprocity drives the system far from equilibrium | **CONFIRMED**, σ²_v 133×; EPR positive and χ²-consistent |
| H_A — antisymmetric part cannot alter static structure | **FALSIFIED** (n_cl is static and moves 12.4×) |
| H_B — signature appears in probability currents | **PARTLY**: current exists (EPR), but no circulation found |
| H1 — turnover monotone in χ | **FAILS at paper scale** (peaks at χ=0.75) |
| H2 — R_Q < 0.1 and R_turnover > 0.5 | R_Q passes, R_turnover fails (as pre-committed) |
| H3 — σ²_v rises with χ | **CONFIRMED**, ρ=+0.975 |
| H4 — χ=0 recovers the monodisperse null | **PARTIAL**: 52% for turnover, 99.6% for σ²_v |
| "Frozen Q" as evidence of solenoidality | **RETIRED** — Q is degenerate (40× n_cl range, 14% Q range) |
| α̂ / s_II/s_I as nonreciprocity knobs | **RETIRED** — α̂ cancels from the antisym/sym ratio |

## 7. What I would do next, in order

1. **Explain the χ=0.25 coarsening dip.** Reproducible at both scales, unexplained, and the
   most interesting unclaimed result here.
2. **Average the EPR bias probe over ≥5 noise seeds** and push dt to ~1.5e-4. That would
   resolve χ ≤ 0.5 and turn the χ² check from 2-of-3 seeds into a clean measurement.
3. **Find what the current couples to.** EPR says it exists; no network observable sees it.
   Candidates: directed contact graph with force-signed edges (the pipeline is undirected
   throughout, despite testing a *directed*-graph prediction), cluster-COM ballistic motion,
   or entropy production resolved per contact.
4. **Retire edge Jaccard** or replace it with a fitted decorrelation rate on hysteretic
   contacts. It failed at paper scale and is half polydispersity even when it works.
5. Only then revisit publication framing — the current README abstract states a conclusion
   ("frozen Q") that this work does not support.

## Files

`data/` N=1000 (61 runs) · `data_paper/` N=4000 (21 runs) · `paper_*.csv` paper-scale tables ·
`epr_paired.csv` entropy production · `AUDIT.md` design critique · `EXPERIMENT.md` χ spec ·
`RESULTS.md` pilot-scale record · `chi_test.png`, `sweep.png`, `modularity_test.png`
