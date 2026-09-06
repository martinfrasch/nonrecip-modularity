# Item 3: raw dyadic data request — scientific case and specific ask

Prepared 2026-09-03. **Updated 2026-09-05: the §7 fallback has been executed** — the analysis is
now validated *and already run* on four published sociomatrices (`DOMINANCE_HODGE_RESULT.md`).
That result does not remove the need for this request; it sharpens it, and changes what we should
ask for. See §4a and the revised §5–§7. Companion analysis: `hodge_dominance.py`.

---

## 1. Why this is the measurement that matters

The main manuscript turns on a distinction it can only test indirectly. A widely held intuition
holds that a nonreciprocal coupling splits into a symmetric sector that selects structure and an
antisymmetric sector that is *solenoidal* — that circulates without a potential. We show this
identification fails in a colloidal system, and note in Section 1 that it was never justified
formally: antisymmetry under particle exchange, violation of Newton's third law,
divergence-freedom of a vector field, and circulation of a probability current are four distinct
properties, and a reciprocity decomposition of a pair force is not a Helmholtz–Hodge
decomposition of a current.

**In a social dominance network that decomposition can be computed exactly.** For a group of n
animals, the net directed interaction flow f_ij = −f_ji is an antisymmetric edge function on a
complete graph. Filling the triangles leaves no harmonic component, so the discrete Hodge
decomposition is exact and two-part:

```
f  =  grad(s)          +  curl
      TRANSITIVE          INTRANSITIVE
      derives from a      no potential exists;
      rank potential s    rock-paper-scissors structure
```

The gradient part is a hierarchy that *is* derivable from a scalar potential — the social analogue
of a conservative force. The cyclic part is genuinely non-gradient. **The cyclic fraction
‖curl‖/‖f‖ is therefore a direct, exact, finite-dimensional measurement of the question our
colloid work can only probe through its consequences: is this antisymmetric coupling a gradient,
or does it circulate?**

This is also the specific measurement that would connect to Section 3.11. Our Vicsek control
showed that coarse-grained circulation detects a *cyclic collective mode*, and that such a mode
can be imposed rather than emergent. A cyclic dominance component is a cyclic collective mode in
a biological social system that nobody imposed.

## 2. Why the published data cannot support it

The open behaviour table (848 rows × 67 columns, MIT licence) is aggregated **per mouse per day**.
It gives each animal's rate of chasing and of being chased, and we have used it productively:
the per-group asymmetries sum to zero to machine precision, and the antisymmetric combination is a
more reproducible individual trait than either direction alone (ICC 0.58–0.78 against 0.18–0.68;
see `FORKOSH_RESULT.md`).

But it does not record **who interacted with whom**. Without ordered pair counts the adjacency
f_ij does not exist, and neither the Hodge decomposition, nor transfer entropy between named
individuals, nor any irreversibility estimate is possible.

## 3. The specific ask, in order of preference

**(a) Minimal and probably sufficient — ordered dyadic counts.** One row per ordered pair per
group per day:

| GroupNumber | Day | actor MouseNumber | target MouseNumber | behaviour | count |
|---|---|---|---|---|---|

for the behaviours already in the aggregate table (aggressive chase, escape, follow, approach,
contact). For 53 groups × 4 days × 12 ordered pairs × ~5 behaviours this is of order 10⁴ rows — a
small file, and derivable from whatever pipeline produced the aggregate columns, since those
columns are row-sums of exactly this table.

**Also needed, and now the decisive quantity:** the per-ordered-pair total interaction count. The
cyclic fraction has a sampling floor set by that count (§4a), so the statistic is uninterpretable
without it — and the published matrices fail precisely here, not on structure. Target density is
≳100 interactions per ordered pair for a marginal test and ≳500 for a decisive one; four days of
continuous video per group should be far above both, which is the whole reason for asking.

**(b) Better — timestamped dyadic events.** One row per event: `(group, day, time, actor, target,
behaviour)`. This additionally permits the sequence to be tested for time-reversibility, an
entropy-production estimate on the dominance-configuration state space, and transfer entropy
between individuals — the same coarse-graining test we apply to colloids and cilia, at a third and
very different scale.

**(c) Fullest — raw trajectories.** `(group, day, frame, MouseNumber, x, y)`. Only needed if
interaction definitions are to be varied; (b) suffices for everything proposed here.

## 4. The analysis is built and validated, not proposed

`hodge_dominance.py` implements the decomposition and validates it on constructed cases:

| case | cyclic fraction | expected |
|---|---:|---:|
| perfectly transitive (ranks recovered exactly) | 0.0000 | 0 |
| pure 4-cycle (ranks flat, as they must be) | 1.0000 | 1 |
| mixtures at cyclic weight 0.25 / 0.50 / 0.75 | 0.316 / 0.707 / 0.949 | monotonic |

It also establishes a methodological point that materially affects the result. The decomposition
must be applied to **log-odds**, not to raw net counts. A Bradley–Terry hierarchy has
logit(p_ij) = r_i − r_j exactly, so its log-odds flow is exactly a gradient; the net-count flow is
a *logistic* function of the rank difference, and a linear decomposition of counts therefore
retains a spurious cyclic fraction that does not vanish with more data:

| interactions per pair | cyclic fraction, raw counts | cyclic fraction, log-odds |
|---:|---:|---:|
| 5 | 0.262 ± 0.093 | 0.288 ± 0.114 |
| 20 | 0.198 ± 0.054 | 0.197 ± 0.079 |
| 100 | **0.157 ± 0.030** | 0.094 ± 0.042 |
| 500 | **0.158 ± 0.013** | **0.044 ± 0.020** |

On raw counts a perfectly transitive group reads as 16% intransitive however much data is
collected. On log-odds the floor falls as it should. Any published claim of intransitive dominance
based on a linear decomposition of counts should be checked against this.

## 4a. The fallback has been run, and it identifies the binding constraint

The `compete` R package (Curley lab) bundles four complete pairwise win–loss sociomatrices. We
decomposed all four (`DOMINANCE_HODGE_RESULT.md`), computing for each a *matched* transitive
floor: fit the rank potential, regenerate a perfectly transitive Bradley–Terry truth at that
matrix's own per-pair counts, re-decompose, and take the resulting distribution as the null.

| dataset | n | interactions | interactions per ordered pair | cyclic fraction | matched floor | z |
|---|---:|---:|---:|---:|---:|---:|
| mouse | 12 | 234 | **3.5** | 0.495 | 0.580 ± 0.054 | −1.6 |
| caribou | 20 | 823 | **4.3** | 0.607 | 0.579 ± 0.028 | +1.0 |
| bonobos | 6 | 739 | 49.3 | 0.341 | 0.228 ± 0.045 | +2.5 |
| people | 6 | 474 | 31.6 | 0.664 | 0.290 ± 0.074 | +5.1 |

Two things follow, and together they are the reason to make this request rather than to treat the
question as settled.

**First, the measurement works.** The two dense matrices have floors near 0.25, resolve cyclicity,
and disagree with each other — so the statistic has real dynamic range when the counts support it.

**Second, the scarce resource is not pairwise structure but interaction density.** The two sparse
matrices sit at 3.5 and 4.3 interactions per ordered pair, where the matched floor is ≈0.58. At
that density the estimator cannot distinguish a perfectly transitive hierarchy from a substantially
cyclic one *in either direction*. The mouse z = −1.6 must therefore be read as "consistent with
transitive, with low power", not as a positive finding of gradient structure.

This is the constraint a round-robin tube test cannot escape: its interaction count per dyad is
fixed by the protocol, typically at single digits, and no amount of additional published matrices
of that design will lower the floor. Only **continuous observation of freely interacting animals**
produces the density required — which is precisely what the Forkosh protocol (four days of
continuous video per group) generates and what the aggregate table's row-sums imply has already
been recorded.

**How much density is needed.** Under a Bradley–Terry truth the log-odds floor falls with the
per-pair count (n = 4 demo, `hodge_dominance.py`): 0.288 ± 0.114 at 5 per pair, 0.197 ± 0.079 at
20, 0.094 ± 0.042 at 100, 0.044 ± 0.020 at 500. A cyclic fraction of 0.15 — a modest departure from
transitivity — is therefore indistinguishable from the floor below ~100 interactions per ordered
pair and cleanly resolved at ~500. That is the number the request should be framed around.

## 5. What we would predict, stated before seeing the data

- **Cyclic fraction small but nonzero, and — critically — resolvable.** Rodent hierarchies are
  usually reported as largely transitive, so we expect the gradient part to dominate. The public
  mouse matrix is consistent with that but cannot demonstrate it (floor 0.580 at 3.5 interactions
  per pair). With continuous-observation counts the floor should fall to ≈0.05–0.10, at which point
  a cyclic fraction indistinguishable from it would be a genuine negative result: social dominance,
  though manifestly nonreciprocal, would be a *gradient* structure — reinforcing the manuscript's
  central point that antisymmetry does not imply circulation. That statement cannot currently be
  made from any published matrix we have found.
- **Cyclic fraction anti-correlated with hierarchy stability.** We found individual asymmetry to be
  a stable trait (ICC 0.58–0.78) while hierarchy steepness is labile (0.24–0.37). If the cyclic
  component is what fluctuates, groups with higher cyclic fraction should show less stable
  steepness across days. This is testable with (a) alone, and the 53-group design gives it power
  that no single published matrix has.
- **Under (b): whether social irreversibility survives coarse-graining.** Colloids: no. Cilia: yes,
  but a synthetic flock can be made to do the same. A social system would be a third point, and the
  first where any cyclic mode present is neither imposed nor a single-particle property.

## 6. Draft of the approach

> Dear Dr Forkosh, Dr Karamihalev and Prof. Chen,
>
> We have been using your published behaviour table (via the IdentityDomains repository) in a study
> of nonreciprocal interactions, and have found something in it we thought you might like: across
> all four directed interaction types, the *asymmetry* between a behaviour and its reverse
> (chase minus escape, follow minus being-followed) is a substantially more reproducible individual
> trait across days than either direction alone — intraclass correlation 0.58–0.78 against
> 0.18–0.68. The per-group asymmetries also sum to zero to machine precision, which is a nice
> internal validation of the rate columns.
>
> We are writing to ask whether ordered pairwise interaction counts — who chased whom, by group and
> day — could be made available. The aggregate columns appear to be row-sums of exactly such a
> table.
>
> Our interest is a Hodge decomposition of the dominance flow, which separates it exactly into a
> transitive component derivable from a rank potential and an intransitive cyclic residual. That
> distinction is central to a question in nonequilibrium physics about whether antisymmetric
> couplings are gradient or circulating, and a small social group is one of the few systems where it
> can be computed exactly rather than inferred.
>
> We have already implemented the analysis and run it on the four complete sociomatrices bundled
> with the `compete` package, and the outcome is what motivates this request. The statistic has real
> range where the data are dense: the two matrices with ~30–50 interactions per ordered pair give
> cyclic fractions well separated from a matched transitive null. But the two matrices built from
> round-robin testing have only 3.5 and 4.3 interactions per ordered pair, and at that density the
> null floor sits at ≈0.58 — high enough that a perfectly transitive hierarchy and a substantially
> cyclic one are indistinguishable. The published rodent matrix is in that regime, so the question
> simply cannot be answered from it.
>
> What breaks that limit is not more matrices of the same design but many interactions per dyad,
> which only continuous observation provides. Four days of continuous video per group is exactly the
> regime where the floor drops far enough for the measurement to mean something, and your 53-group
> design would additionally let us test whether cyclicity tracks the stability of the hierarchy
> across days — which the single published groups cannot address at all.
>
> We would of course propose collaboration and co-authorship on anything arising, rather than merely
> a data transfer, and we are happy to send the current analysis and the four-dataset result first
> if that is useful.

## 7. Status of the fallback

**Executed 2026-09-03/05.** The literature/package search returned four complete pairwise
sociomatrices, all four have been decomposed, and the result is in §4a and
`DOMINANCE_HODGE_RESULT.md`. The fallback therefore no longer functions as a substitute for this
request — it *is* the argument for it, because it quantifies exactly why the published data cannot
answer the question and what property new data would need.

Remaining fallback options, in order of value:

1. **Other continuous-observation rodent datasets.** Curley-type protocols (1–3 h daily for 21
   days) reach interaction densities comparable to the bonobo/people matrices. Any such dataset with
   ordered dyadic counts preserved would support the measurement without a new request.
2. **Dense non-rodent social data.** The bonobo and people matrices show the estimator resolving
   cyclicity; more groups at that density, in any species, would establish whether the transitive
   result generalises.
3. **Aggregation of many sparse matrices from the same protocol.** Pooling does not lower the
   per-group floor, but a meta-analysis of per-group cyclic fractions against their own matched
   floors would have more power than any single group. This is the cheapest remaining option and
   requires no new data.

Caveats on everything above are in `DOMINANCE_HODGE_RESULT.md`: four single unreplicated groups,
one per dataset; the floor is computed under a Bradley–Terry null; and the cyclic fraction is a
structural statistic, not a thermodynamic one — it answers "is this antisymmetric coupling a
gradient?", not "how much does this system dissipate?".
