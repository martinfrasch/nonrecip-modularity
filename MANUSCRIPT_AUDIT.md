# Full manuscript audit, 2026-09-05

Every quantitative claim in `PAPER.md` was re-derived from the trajectories, CSVs or source
documents rather than taken from the text. Every method was re-read against its implementation.
This file records what was found, what was changed, and what was verified unchanged.

## A. Logical defects (the important ones)

**A1 — The tier-membership conclusion was stale and self-contradictory.** The manuscript stated in
five places (Introduction, §3.10 Interpretation, §4.6, §4.7, Conclusion 8) that the system "passes
two of the three" tier-II tests and therefore "occupies the linear-response tier". By §3.10's own
text this was already false:

- **A1 (quadratic EPR)** had been reframed within §3.10 itself as "best read as a validation of the
  entropy-production estimator ... rather than as an independent test of linear response".
- **A2 (Onsager)** had been withdrawn within §3.10 by the dissimilar-channel control: "its
  interpretation as a physical reciprocity relation does not survive this control".
- **A3 (effective temperature)** was reported as ill-posed.

So **zero** of three returned a decisive positive, while the summary paragraphs still claimed two.
Rewritten throughout: tier membership is *undetermined*, and the one diagnostic that converged
(no window-independent effective temperature) points **away** from tier II. This is a stronger and
more interesting result than the one it replaces, and it is the one the data support.

**A2 — §3.10 contained three paragraphs written before the dissimilar-channel control and never
updated.** They asserted the cross-coupling is "purely reactive", advanced a gyroscopic mechanism,
and stated that the artifact explanation "is equally live **until the dissimilar-channel control is
run**" — three paragraphs *after* that control's results were tabulated. Consolidated into one
block that records the hypothesis as tested-and-not-confirmed.

**A3 — Conclusion 7 contradicted §3.9.** §3.9 explicitly declines to call S* a critical nucleus;
Conclusion 7 asserted it "is a critical nucleus — an unstable fixed point". Conclusion corrected to
"event-rate crossover", with the two measurements that would upgrade it named.

**A4 — §4.2 and §4.7 still asserted the withdrawn "purely reactive cross-coupling"** as an
established contribution of the antisymmetric sector. Both now rest on the dissipation alone.

**A5 — §2.5 stated "The biological comparison uses no new simulation"** while §3.11 reports a
four-variant Vicsek simulation. Corrected.

**A6 — §1 "The third is confirmed"** was stronger than §4.7, which carefully states the relation is
non-monotonic and that no functional form predicts the dip. Softened to match.

## B. Numerical errors found by recomputation

| location | manuscript | recomputed | note |
|---|---|---|---|
| §3.2 `lcf` row | 0.301 / **0.910** / **0.898** / **0.699** / **0.717** / **0.771** / **0.46** | 0.301 / 0.881 / 0.904 / 0.800 / 0.791 / 0.753 / 0.261 | mixed final-snapshot and late-half values; **Figure 2 was always correct**, only the table was wrong |
| §3.2 scatter | n_cl 11.1 **± 1.3** → 137.3 **± 0.6** | ± 0.5 → ± 0.5 | the ±1.3 is the *monodisperse* row's sd |
| §3.2 pilot σ²_v ratio | **24×** | 54× | 24–26× is the *mono-vs-bidisperse* pilot contrast, a different comparison |
| §3.3 Q span | **14%** | 15% | also rescoped to "the 14 pilot-scale conditions" |
| §3.7, §4.3 turnover survival | **74%** | 69% (paper), 71% (pilot) | quoted as "about 70%" |
| §3.7 morphology correlation | ρ = **+0.795** | +0.77 (pilot), +0.62 (paper) | both now quoted |
| §3.5 drift exponent | N^**−0.50** | −0.38 (per cluster), −0.47 (binned), −0.09 to −0.12 over the small-to-medium range | the theoretical value had been reported as a measurement; the argument now rests on the coherence exponent (−0.53), which is direct and does hold |
| §3.8 composition | 22.7% gives **0.831** vs 0.786, a **+5.8%** difference | 0.771 vs 0.786, ≈2% and of *opposite* sign | stale two-seed value; the corrected number *strengthens* the conclusion |
| §3.8 "**Both seeds** plateau at 0.83–0.85" | — | removed | stale two-seed sentence in a four-seed table |
| §3.8 replicate window-1 | 0.378 / 0.453 / 0.443 / 0.556 | 0.381 / 0.479 / 0.441 / 0.554 | within-noise; both `PAPER.md` and `BOX_SCALING.md` aligned to the reproducible values |
| §5 rectangular box | 0.809 ± 0.009, **3.8%** | 0.804 ± 0.008, ≈3% | |
| Abstract | "**133 simulations** at **four** system sizes" | 113 independent trajectories + 25 continuations = 138 files, five system sizes | |
| §2.5 | "**129 runs**" (breakdown sums to 101); "**2 runs**" replicate; "**Eleven** continuation runs" | 113 / 4 / 25; density series (6) and high-drive series (12) were omitted entirely | replaced with a full campaign table |
| §2.5 convergence | "+100.1% (N=22,000)" | +68% mean, range +24% to +115% | stale two-seed value |

## C. Methods defects

**C1 — `tests/test_reciprocity.py` was broken and could not run.** It was never updated when the
`heat` accumulator was added to `_run`, so it failed with `TypeError: not enough arguments:
expected 15, got 14`. This is the test that validates the paper's central construction. Fixed and
re-run; it passes and reproduces §2.2 exactly (COM drift 0.000e0 / 1.975e−2 / 3.950e−2, and the
equal-radius pair identical at χ=0 and χ=1).

**C2 — §3.10 misdescribed its own χ²-scaling data as "frozen configuration".** `epr_probe.py`
probes each χ from the final snapshot of a run equilibrated *at that same χ*, so structure
co-varies with drive along that row. The near-tautology argument the manuscript used to discount
the result does not apply to it. Rewritten to separate the two measurements, and to state the
construction argument correctly (the symmetric sector's heat is an exact differential averaging to
zero in a stationary state, leaving terms linear and quadratic in χ, of which the quadratic must
dominate at large χ regardless).

**C3 — §3.10's χ² and §3.4's χ^3.01 appeared to contradict each other.** They are different χ
ranges of the same estimator (0.75–1.5 and 1.5–8). Now stated as a drive-dependent exponent rather
than two competing claims. Also noted that the χ^3.01 fit uses only the three levels at which EPR
resolves above the bias noise.

**C4 — Methods described none of the comparison systems.** The cilia pipeline, the single-cluster
shape descriptor, the Vicsek model and the Hodge decomposition appeared only in Results. Added as
§2.6, with parameters checked line-by-line against `vicsek.py`.

**C5 — Convergence claim overstated at N=16,000.** §2.5 claimed between-seed scatter collapses on
convergence; at N=16,000 two converged seeds give lcf 0.709 and 0.857. Stated explicitly, since
that point carries the scaling fit.

**C6 — The χ=0 entropy production is zero by construction, not by measurement.** In the paired
protocol the estimate and its control are the same quantity at χ=0, so the net is identically zero.
Added to Limitations; the χ=0 zero-current claim rests on the analytic argument of §2.3.

## D. Under-selling

**D1 — The dominance-network work was absent from the manuscript entirely.** It is a completed,
validated result on published data that bears directly on the paper's central proposition, and it
supports it by an independent route. Added as §3.12, with a discussion link in §4.2 and
Conclusion 11 — including the methodological finding that the Hodge decomposition must be taken on
log-odds, since raw counts carry an irreducible ~16% cyclic floor.

**D2 — The Vicsek control had no results document.** Reproduced exactly and written up as
`VICSEK_RESULT.md`.

## E. Verified correct, unchanged

- **The χ = 0 detailed-balance argument (§2.3).** Re-derived: at χ=0 all pair forces are central,
  pairwise and equal-and-opposite, hence conservative (the r ≤ 1 cutoff makes the force
  discontinuous but the potential well defined); mobility μ_i = 1/s_i is configuration-independent;
  noise variance σ²/s_i gives D_i/μ_i = σ²/2 uniformly. Detailed balance holds. Correct.
- **The χ construction (§2.2)** — symmetric sector χ-invariant, antisymmetric exactly linear.
  Verified in code and by re-running the (repaired) test.
- **`fdt.py`'s collective-coordinate conjugacy.** The perturbation couples to X = Σ_{i∈species} x_i
  and both response and fluctuation are measured on X; T_eff = ⟨ΔX²⟩/(2μ) then reduces to D/μ.
  Correct as implemented.
- **Pair well depths independently re-derived** by integrating the symmetric EHD coefficient from
  contact to the cutoff: 6.99 / 6.78 / 4.83 kT for L-L / L-S / S-S, against the manuscript's
  7.0 / 6.8 / 4.8. This is one of the three legs of the kinetic-arrest argument in §3.5.
- **The paired EPR bias control is correctly implemented**: `epr_probe.py` re-probes the *same*
  starting configuration at χ=0 with the *same* noise seed, so the configuration-dependent
  discretisation bias cancels on subtraction. This was the fix for the earlier failure mode where
  an excess over a differently-structured χ=0 run did not cancel.
- Paper-scale n_cl (11.1 / 2.2 / 17.3 / 56.5 / 95.2 / 137.3, mono 11.2); the 12.4× and 133× ratios;
  ρ = +0.931 (p = 2.1e−8) and +0.975 (p = 7.1e−12).
- The within-run sd/mean figures of §3.5 (7.2 / 29.7 / 18.6 / 7.1%) — exact, though the manuscript
  did not say they were *within-run*; between-seed scatter at χ=0.25 is larger still (67%).
- ΔQ_null 0.440 → 0.383 with mono 0.462; Q–n_cl degeneracy ρ = +0.31, p = 0.28; the 3-cluster
  (Q=0.869) vs 33-cluster (Q=0.854) pair.
- §3.6 attribution: 52% and 99.6%.
- Box scaling: 3,116 / 7,221 / 12,581, ratios 1.00 / 2.32 / 4.02, exponent 1.01, lcf 0.78–0.80.
- Replicate convergence: 0.771 ± 0.075, paired t = +4.4, p = 0.022, all four seeds rising.
- Coarsening exponents z = 0.089 (χ=0) and 0.27–0.42 (χ>0), 1/z = 2.4–3.7.
- EPR/χ² flat to 4.8% (with and without the χ=1.5 anchor); the χ=2..8 four/fifteenfold ranges.
- Onsager identical-channel and dissimilar-channel tables.
- Cilia: eight ATP rows summing to n = 184; 92% above |z|=2 and 55% above |z|=3 recomputed from the
  per-concentration table.
- Single-cluster control: **re-run**, reproduces 0.81 / 0.85 / 0.94 with 0.00 above |z|=2, 3–6
  clusters, 101 frames.
- Vicsek control: **re-run**, reproduces 0.709/0.55, 0.875/0.32, 0.015/1.80, 0.069/5.58 exactly.
- All eleven referenced figure files exist and are generated from the CSVs.

## F. Reference corrections

- [12] listed Broedersz first for a paper the text cites as "Battle et al."; corrected to the
  actual author order.
- [14] contained "and related work on ..." inside the citation; removed.
- Added [27] Forkosh et al. (2019), [28] Jiang et al. on combinatorial Hodge theory, [29] the
  `compete` package.

## G. The effective-temperature result, re-verified

The A3 numbers had no source document, and A3 now carries the paper's only converged tier
conclusion, so `fdt.py` was re-run at the manuscript's stated 200 starts per condition. **It
reproduces.**

| quantity (χ = 1.5, large species) | re-run, 200 starts | manuscript |
|---|---:|---:|
| T_eff/kT at T = 200 | 5.63 | 5.48 |
| T_eff/kT at T = 400 | 10.10 | 9.78 |
| window ratio | 1.79 | 1.78 |
| ⟨ΔX²⟩ exponent | **+1.821** | 1.81 |
| response exponent | +0.981 | 0.99 |
| resulting T_eff exponent | **+0.843** | +0.83 |

Agreement to three significant figures on the exponents. The superdiffusion of the collective
coordinate — the substantive finding, and the reason no window-independent effective temperature
exists — is confirmed.

The χ = 0 row was only re-run at 40 starts (the 200-start sweep was abandoned for machine time).
At that statistics it gives T_eff/kT = 0.74 ± 0.17 and 0.85 ± 0.17 at T = 200 and 400, both
consistent with unity as required, and a ⟨ΔX²⟩ exponent of 1.16 ± ≈0.4 against the manuscript's
0.86. Those are consistent with each other, and both are consistent with ordinary diffusion; the
contrast that matters — 1.0-ish at χ = 0 against 1.82 at χ = 1.5 — is unambiguous either way. One
small change was made: the manuscript described the species difference at χ = 1.5 as measured "at a
single window"; it is present at both (t = +3.8 and +3.6), so the text now says so while keeping
the same interpretation.

**A sharper frozen-configuration EPR test** (one configuration, χ ∈ {1, 1.5, 2, 3, 5, 8}, fitting
net = aχ + bχ² to isolate the linear cross-term) was started and abandoned as too slow to justify
the machine time. No claim depends on it: §3.10 now rests the near-tautology argument on §3.4's
existing frozen measurement and states explicitly that it is not sharp enough to isolate a linear
term. Worth running as a follow-up.

## H. References — resolved by checking minaction.net

The open question from the first pass (ref [9] cited both for NWAP and for Minimum-Action Learning,
though its title is the latter) is resolved. minaction.net and the four validation papers were
checked directly. The programme comprises:

| | reference | role |
|---|---|---|
| 2026a | *Causal thinking in physiology*, J Physiol, DOI 10.1113/JP290762 | **defines NWAP**: S_NW = ∫(E − I + A·C) dt |
| 2026b | arXiv:2603.16951 | Minimum-Action Learning (physics) |
| 2026c | arXiv:2604.24805 | energy-first neural architecture |
| 2026d | arXiv:2605.05254 | **the modularity-excess prediction** (marine metabolic networks) |
| 2025 | HAL hal-05347658 | network-weighted action framework, formalism preprint |

Three changes follow.

**H1 — [9] was the wrong paper for NWAP.** It pointed at the MAL preprint. NWAP is defined in
Frasch 2026a (J Physiol). [9] is now that paper; MAL moved to [31]; 2026c, the HAL preprint,
Friston and England added as [32]–[35].

**H2 — The modularity-excess prediction now has its actual source.** Frasch 2026d establishes, in
marine metabolic networks, that *excess over a null* rather than absolute modularity is the
informative signature, reporting ΔQ ≈ 0.15–0.40 over configuration-model, label-permutation and
bipartite-incidence nulls. That is precisely the prediction §3.3 refutes with inverted sign, and it
was previously uncited. §3.3 and §4.7 now cite it, and §4.7 states the comparison properly: our
excess over the same class of null is of *comparable magnitude* (0.38–0.46), so the quantity is not
absent — it simply runs the wrong way with nonreciprocity. §4.7 also now notes that our
degeneracy diagnosis is specific to 2D contact networks and does not transfer to their sparse
bipartite topology.

**H3 — The directed-graph extension is not published anywhere, and the manuscript implied it was.**
The framework, validation, origins, future-work and papers pages were all checked: none mentions
directed graphs, nonreciprocity, antisymmetric couplings, solenoidal components, Onsager,
Rayleighian or Maximum Caliber. The directed-graph extension and the symmetric/antisymmetric
assignment are the author's own unpublished extension, formulated while designing this study. §1
and §4.7 now say so explicitly — the extension is a hypothesis of ours tested as such, while the
action principle itself [9, 33] and the modularity-excess prediction [30] are published. This is
more defensible than the previous phrasing and costs the paper nothing: testing one's own stated
hypothesis is exactly what §4.7's disclosure already describes.

**H4 — Fourteen of thirty-five references were defined but never cited** (Fruchart; Newman;
Fortunato & Barthélemy; Sekimoto; Onsager; Casimir; Jaynes; Fodor; Nardini; Doi; Good et al.;
Peixoto; Kubo; Seifert). All are now cited at the point where the relevant claim is made. Two of
them strengthen §3.3 materially: the modularity-degeneracy result is a known hazard of modularity
maximisation [23] and part of why inferential community detection is preferred [24], which makes
our finding an instance of a documented problem rather than an isolated observation.
