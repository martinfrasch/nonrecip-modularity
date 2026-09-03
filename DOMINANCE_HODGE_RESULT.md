# Hodge decomposition of published dominance matrices — first pass, no data request required

Item 3 fallback, run 2026-09-03. It turns out the fallback is available immediately: the
`compete` R package (Curley lab, [github.com/jalapic/compete](https://github.com/jalapic/compete))
bundles four **complete pairwise win–loss sociomatrices**, which is exactly the input the Hodge
analysis needs. Exported to `data_dominance/`.

## Result

Cyclic fraction of the log-odds flow, against a sampling floor computed by fitting the rank
potential, regenerating a *perfectly transitive* Bradley–Terry truth at that matrix's own
per-pair interaction counts, and re-decomposing:

| dataset | n | interactions | cyclic fraction | transitive floor | z |
|---|---:|---:|---:|---:|---:|
| mouse | 12 | 234 | 0.495 | 0.580 ± 0.054 | **−1.6** |
| bonobos | 6 | 739 | 0.341 | 0.228 ± 0.045 | **+2.5** |
| caribou | 20 | 823 | 0.607 | 0.579 ± 0.028 | +1.0 |
| people | 6 | 474 | 0.664 | 0.290 ± 0.074 | **+5.1** |

**The mouse hierarchy is, if anything, *more* transitive than a matched Bradley–Terry model
(z = −1.6).** Its antisymmetric social coupling is a gradient: derivable from a scalar rank
potential, with no detectable circulating component. That is consistent with these CD1 groups
being described as forming "extremely linear" hierarchies, and it is the same conclusion the
colloid work reaches by a different route — **antisymmetry does not imply circulation.**

The human dataset is strongly cyclic (z = +5.1) and the bonobos modestly so (z = +2.5). We do not
interpret the species pattern: these are single groups, unreplicated, and we do not know what
competition the human data records.

## The binding constraint is interaction density, not data availability

| dataset | pairs | interactions per pair | floor |
|---|---:|---:|---:|
| mouse | 66 | **3.5** | 0.580 |
| caribou | 190 | **4.3** | 0.579 |
| bonobos | 15 | 49.3 | 0.228 |
| people | 15 | 31.6 | 0.290 |

The two sparse matrices have floors near 0.58 — **they cannot resolve cyclicity at all**, and the
mouse z = −1.6 should be read as "consistent with transitive, with low power" rather than as a
positive finding. The two dense matrices have floors near 0.25 and can resolve it, and they
disagree with each other, which shows the measurement has range when the counts support it.

This reframes the request in `REQUEST_RAW_TRACKING.md`. The scarce resource is not pairwise
structure — that is published — but **many interactions per pair**, which only continuous
observation provides. Forkosh et al. record four days of continuous video per group; Curley-type
protocols observe one to three hours daily for 21 days. Either yields orders of magnitude more
interactions per dyad than a round-robin tube test, which is precisely what is needed to push the
floor down.

The ask should therefore emphasise **interaction counts per ordered pair**, and note that the
analysis is already validated and already run on sparser public data, with the density limit
quantified above.

## Caveats

- Four single groups, one per dataset; no replication within any species.
- The floor is computed under a Bradley–Terry null. A different transitive generative model would
  shift it, though the log-odds scale removes the main nonlinearity confound (see
  `hodge_dominance.py`).
- Nothing here measures dissipation or irreversibility. The cyclic fraction answers "is this
  antisymmetric coupling a gradient?", which is a structural question, not a thermodynamic one.
