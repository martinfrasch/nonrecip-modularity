# Item 3: raw dyadic data request — scientific case and specific ask

Prepared 2026-09-03. Companion analysis: `hodge_dominance.py` (implemented and validated).

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

**Also needed:** the per-pair total interaction count. The validation below shows the cyclic
fraction has a sampling floor that depends on it, so the number is not interpretable without it.

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

## 5. What we would predict, stated before seeing the data

- **Cyclic fraction small but nonzero.** Rodent hierarchies are usually reported as largely
  transitive, so we expect the gradient part to dominate. A cyclic fraction indistinguishable from
  the matched sampling floor would be a clean negative result and would say that social dominance,
  despite being manifestly nonreciprocal, is a *gradient* structure — reinforcing the manuscript's
  central point that antisymmetry does not imply circulation.
- **Cyclic fraction anti-correlated with hierarchy stability.** We found individual asymmetry to be
  a stable trait (ICC 0.58–0.78) while hierarchy steepness is labile (0.24–0.37). If the cyclic
  component is what fluctuates, groups with higher cyclic fraction should show less stable
  steepness across days. This is testable with (a) alone.
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
> We would like to ask whether ordered pairwise interaction counts — who chased whom, by group and
> day — could be made available. The aggregate columns appear to be row-sums of exactly such a
> table. Our specific interest is a Hodge decomposition of the dominance flow, which separates it
> exactly into a transitive component derivable from a rank potential and an intransitive cyclic
> residual. That distinction is central to a question in nonequilibrium physics about whether
> antisymmetric couplings are gradient or circulating, and a small social group is one of the few
> systems where it can be computed exactly rather than inferred.
>
> The analysis is implemented and validated, and we would of course propose collaboration and
> co-authorship on anything arising, rather than merely a data transfer. We are happy to send the
> current analysis and results first if that is useful.

## 7. Fallback if the data are unavailable

Several published rodent-dominance datasets report full pairwise win/loss matrices (tube tests,
warm-spot competition, home-cage observation). These are smaller and lack the longitudinal
behavioural breadth, but they support the Hodge measurement directly and would establish the
cyclic-fraction result independently of any single group's data. That is the route to take if this
request is declined or slow, and it is worth starting a literature search for such matrices in
parallel rather than sequentially.
