# Design audit: do the tests actually test the hypotheses?

Written 2026-08-22, before the paper-scale repeat. Every number below comes from the 61
runs already on disk; nothing here needs new simulation.

## The hypotheses, as stated in README.md

> **H_A** "the antisymmetric (nonreciprocal) part of the coupling is solenoidal — it cannot
> alter a static structural observable"
> **H_B** "its entire signature appears in probability currents: frozen Q, circulating partition"
> **H_C** "turnover rate at fixed Q tracks the antisymmetric/symmetric coupling ratio"

**Summary of the audit: H_A is falsified by data already in hand; H_B has never been
measured; H_C was tested with a knob that could not vary the ratio (fixed in §4 of
RESULTS.md) and then with one that could, but using a contaminated outcome variable.**

---

## Finding 1 — Q is a saturated observable, so "frozen Q" is close to vacuous

Across all 14 conditions run, cluster count spans **40×** (1 → 40) while Q spans **14%**,
and the two are uncorrelated: Spearman(n_cl, Q) = **+0.31, p = 0.28**.

| state | n_cl | Q |
|---|---:|---:|
| single condensed blob (s_ratio 1.5) | 1.0 | 0.772 |
| coarsened, reciprocal (χ=0) | 3.4 | 0.869 |
| fragmented (χ=1.5) | 33.1 | 0.854 |
| most fragmented (α̂=0.003) | 40.2 | 0.890 |

A 3-cluster state and a 33-cluster state are **indistinguishable in Q**. Newman modularity
of a 2D contact network at this density sits near 0.86 essentially regardless of morphology
— Louvain's resolution limit (~√(2E) ≈ 60 nodes here) merges anything smaller, and the
partition is dominated by subdividing whichever blob is largest.

So Q is not evidence of anything being "frozen". It is an observable with almost no
dynamic range on this system. **Any** physics would have left it near 0.86.

## Finding 2 — H_A is falsified: the antisymmetric coupling *does* move static structure

`n_cl` and `lcf` are computed from a **single snapshot**. They are static structural
observables by construction. Across the χ sweep, with the symmetric coupling held bit-identical:

| χ | 0 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 |
|---|---:|---:|---:|---:|---:|---:|
| n_cl | 3.4 | 2.0 | 3.5 | 12.1 | 23.3 | 33.1 |
| lcf | 0.69 | 0.78 | 0.94 | 0.89 | 0.84 | 0.79 |

χ=0 → 1.5 changes n_cl by **9.8×**, Welch t = −32.2, **p = 1.4e-09**. Turning on the
antisymmetric part, and nothing else, restructures the system by an order of magnitude.

This was visible in the pilot's own table (`n_clusters: ~5 vs ~21`, a 4× difference sitting
two rows below the Q row) and was not acted on. The inference "Q didn't move, therefore no
static observable can move" generalised from one insensitive observable to all of them.

**The premise is also wrong in general.** A non-gradient force field does not preserve the
stationary density; it does so only in the special case ∇·(ρ_eq v) = 0. And the paper this
model comes from reports *arrested coarsening* — a nonreciprocity-induced change in steady-state
structure. H_A contradicts the paper it is built on, and the data side with the paper.

## Finding 3 — no measurement in this pipeline is a current

χ=0 is a **rigorous equilibrium reference**, not merely a control:

- all forces are central, pairwise, equal-and-opposite → conservative
- mobility μ_i = 1/s_i with noise variance σ²/s_i → D_i/μ_i = σ²/2, uniform
- fluctuation–dissipation holds ⇒ detailed balance ⇒ **probability current is exactly zero**

Yet at χ=0:

| observable | χ=0 (zero current) | χ=1 | equilibrium fraction |
|---|---:|---:|---:|
| edge Jaccard turnover | 0.273 | 0.369 | **74%** |
| σ²_v | 1.45e-9 | 3.47e-8 | 4% |

**74% of the headline "turnover" is equilibrium churn carrying no current at all.** Both
metrics are activity measures — magnitude of motion — not currents. A hotter equilibrium
system scores higher on both while remaining perfectly reversible. σ²_v is the cleaner of
the two but is still not a current.

## Finding 4 — the "circulating partition" claim is unsupported

The first actual test of circulation in this project: the signed area rate in observable
planes, which is exactly zero under detailed balance for *any* observable pair and nonzero
iff time-reversal is broken. Measured in four planes — (n_cl,lcf), (Q,n_cl), (lcf,edge_jac),
(n_cl,edge_jac) — at every χ:

**Every value is consistent with zero, and none is distinguishable from the χ=0 equilibrium
null (all p > 0.29).**

This is a null result with weak power (5 seeds × ~50 snapshots), so it does not *disprove*
circulation. It establishes that the project's central metaphor has never been measured, and
that the first attempt to measure it finds nothing. The paper-scale runs (nsnap=200) double
the statistics for this test.

## Finding 5 — edge turnover is largely a structural readout

Across the 13 bidisperse conditions, Spearman(n_cl, edge_jac) = **+0.795, p = 1.2e-03**.
Turnover rises and falls with cluster morphology. It is not an independent dynamical channel;
much of its variation is the same information n_cl already carries.

## Finding 6 — edge Jaccard is not a rate, and has a threshold-flicker floor

Jaccard vs lag, late window:

| χ | lag1 | lag2 | lag4 | lag8 | lag16 | lag32 |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.274 | 0.292 | 0.314 | 0.352 | 0.405 | 0.475 |
| 1 | 0.343 | 0.402 | 0.479 | 0.571 | 0.672 | 0.778 |

The lag-1 value does not start near zero — at equilibrium it is already 0.274. A large part
of it is **hard-cutoff flicker**: particles sitting near r = 1.1(s_i+s_j) whose edge blinks
with thermal jitter. That floor is density- and morphology-dependent, which is a second route
by which structure contaminates the metric. A defensible replacement is a fitted decorrelation
*rate* from the multi-lag curve, with a hysteretic (two-threshold) contact definition.

## Finding 7 — σ²_v is not a velocity variance

As implemented it is `Var(|Δx|/Δt_snap)` with Δt_snap = 1000 nondim units = **20,000
integration steps**, over which a particle moves ~0.6 radius. It is the variance of *speed*
(not of velocity) over a coarse-grained interval — a cage-escape/mobility probe.
⟨|v|²⟩ differs from it by 1.7×. Whatever it measures, it is not the paper's instantaneous
σ²_v, so cross-comparison to Hara et al. Fig. 2d is not on equal terms.

## Finding 8 — H1's significance was overstated

Every χ arm shares its seed, giving identical initial configuration **and** identical noise
stream. The 30 runs are 5 blocks of 6 coupled runs, not 30 independent samples, so the
unpaired Spearman (ρ=+0.860, p=1.1e-09) is anticonservative. The correct blocked test:
within-seed ρ = +1.00, +0.66, +0.77, +0.94, +0.94 — all 5 positive, **sign-test p = 0.06**.

The trend is real and the common-random-numbers design is a good one; the p-value was not.

## Finding 9 — smaller items

- **mono vs bi is confounded by N** (500 vs 1000 particles). The "4× more clusters" is 2× in
  mean cluster size, which is the size-independent statement.
- **ARI noise floor** is estimated from a *single* snapshot but compared against a mean over
  ~50 snapshot pairs, and it scatters 0.35–0.87 across runs. Inconsistent comparison.
- **The analysis graph is undirected** (`nx.Graph`) throughout, while the stated target is a
  "directed-graph prediction". The directed content of NWAP is untested by this pipeline.
- `json` is imported and unused in `analyze.py`.
- `pos % Lbox` can return exactly `Lbox` for tiny negative values, which `cKDTree(boxsize=)`
  rejects. Not triggered in 61 runs; latent.

---

## What the evidence actually supports

Restated so each claim matches what was measured:

1. **Nonreciprocity causes arrested coarsening.** n_cl 3.4 → 33.1 across χ with everything
   else fixed, p = 1.4e-09. Strong, clean, and it reproduces the paper's central result with
   a control the paper does not have.
2. **Nonreciprocity raises activity.** σ²_v ×24 from χ=0 to χ=1, with only a 4% equilibrium
   baseline.
3. **Newman modularity is insensitive to all of it** — but that is a fact about Newman
   modularity on 2D contact networks, not about solenoidality.
4. **Nothing has been shown about probability currents**, in either direction.

## Recommended changes before the next round

| priority | change | why |
|---|---|---|
| high | report n_cl / mean cluster size / lcf as the primary structural outcomes; demote Q | Q has no dynamic range here (F1) |
| high | measure irreversibility directly (EPR or area rate at proper statistics) | H_B is otherwise untestable (F3, F4) |
| high | retire H_A as stated; the honest claim is *Q* is insensitive, not *structure* is | F2 |
| medium | replace edge Jaccard with a fitted decorrelation rate + hysteretic contacts | F6 |
| medium | analyse the χ sweep as the blocked design it is | F8 |
| medium | define σ²_v at the integration timestep, or rename it | F7 |
| low | match N between mono and bi, or report mean cluster size | F9 |
| low | build the directed graph if the directed prediction is to be tested | F9 |

## Status of the entropy-production instrument

A Stratonovich heat accumulator (Σ F_i ∘ dx_i, midpoint rule) was added to the kernel and
**failed validation**: at χ=0, where entropy production is exactly zero, it returns a large
positive constant that is independent of χ. The estimator subtracts two nearly-cancelling
O(μ|F|²) terms and the trapezoidal residual — which scales as √dt and is proportional to run
length, not to physics — swamps the signal at the production timestep.

It is retained in `simulate.py` (it costs ~2% runtime and is correct in the dt→0 limit) but
**must not be used at dt=0.05**. Measuring EPR properly needs its own short, small-dt protocol
started from an equilibrated configuration, not the long production runs. Not yet built.
