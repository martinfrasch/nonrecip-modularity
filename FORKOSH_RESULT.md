# Directed social asymmetry as a stable trait: analysis of Forkosh et al. open data

Item 2 of the external-systems programme (`EXTERNAL_SYSTEMS.md`). Data: Forkosh et al.,
*Nat. Neurosci.* **22**, 2023 (2019), behaviour table from
[github.com/OrenForkosh/IdentityDomains](https://github.com/OrenForkosh/IdentityDomains) (MIT).
212 mice in 53 groups, four days each, 848 rows × 67 readouts. Analysis: `forkosh_asymmetry.py`.

## Why this system at all

Group-housed mice are a nonreciprocal interacting system in the strict sense: if A chases B, B
does not chase A. Dominance is the canonical biological instance of broken action–reaction
symmetry. The published table records **both directions** of several social interactions per
animal per day, which is the antisymmetric sector of a social coupling, measured.

We do **not** claim a connection at the level of the mathematics. The identity-domains framework
maximises a between- over within-individual variance ratio — a Fisher-type generalised eigenvalue
problem. It is a variance decomposition, with no action, path measure, dissipation or
time-reversal operation, and the resemblance to a variational principle is superficial.

## 1. A conservation check that validates the reading

Within a closed group every chase is another animal's escape, so the per-group asymmetries must
sum to zero. They do, exactly:

| interaction | mean \|per-group sum\| | typical individual \|asymmetry\| | ratio |
|---|---:|---:|---:|
| aggression (chase − escape rate) | **0.0000** | 1.107 | 0.000 |
| following (follow − being-followed rate) | **0.0000** | 3.371 | 0.000 |
| chase per contact | 0.0188 | 0.034 | 0.559 |
| follow per contact | 0.0574 | 0.101 | 0.568 |

The two *rate* measures are conserved to machine precision, confirming they are genuinely the two
directions of one interaction within a closed group. The two *per-contact fraction* measures are
not conserved, as expected since each is normalised by that animal's own contact count. This
distinction matters for any subsequent use: the rates carry the conservation structure, the
fractions do not.

## 2. The antisymmetric combination is a better-defined trait than either direction

Intraclass correlation across the four days — the fraction of variance that is between-individual
rather than day-to-day — applied to the asymmetry and, separately, to each direction alone:

| interaction | asymmetry ICC | forward alone | reverse alone |
|---|---:|---:|---:|
| aggression | **0.576** | 0.330 | 0.175 |
| following | **0.619** | 0.378 | 0.214 |
| chase per contact | **0.719** | 0.583 | 0.491 |
| follow per contact | **0.781** | 0.676 | 0.612 |

**In every case the antisymmetric combination is more reproducible than either component.** The
difference between chasing and being chased is a more stable property of an animal than either
chasing or being chased. For aggression the gain is substantial: ICC 0.576 against 0.330 and
0.175.

This is the observation of interest. In the colloid system the antisymmetric sector likewise
proves to be the well-conditioned object: the reciprocity parameter χ isolates a quantity that
controls cluster statistics, whereas the individually-natural parameters (coupling strength,
steric size ratio) do not vary it at all. Finding the same structural pattern — the antisymmetric
combination carrying more signal than its components — in a biological system on aggregated
behavioural data is suggestive of why a symmetric/antisymmetric decomposition is a natural
coordinate choice for nonreciprocal systems generally.

We are careful about the strength of that statement. This is a shared structural pattern, not
evidence for a variational principle. Nothing here measures dissipation, irreversibility or a
current.

## 3. Group-level hierarchy is much less stable than the individual trait

| interaction | mean steepness (within-group sd of asymmetry) | group-level ICC across days |
|---|---:|---:|
| aggression | 1.506 | 0.243 |
| following | 4.580 | 0.366 |

So individual asymmetry is a reasonably stable trait (ICC 0.58–0.78) while the *steepness* of the
hierarchy it produces is not (0.24–0.37). The animals keep their dispositions; the group structure
those dispositions generate is comparatively labile. That asymmetry between levels is itself worth
noting, and is the kind of cross-scale observation the vertical/horizontal-scales programme is
about.

## What cannot be done with these data, and what would be needed

The table is aggregated **per mouse per day**. It gives each animal's rate of chasing and of being
chased but not *who chased whom*, so the directed adjacency G_ij is not recoverable. Consequently:

- no directed-network analysis, no transfer entropy between named individuals;
- no test of whether the dominance structure is transitive or **cyclic** — cyclic dominance being
  the genuine social analogue of circulation, and the measurement that would connect this system
  to Section 3.11 of the main manuscript;
- no irreversibility or entropy-production estimate.

All three require the raw tracking data. Given the conservation check above and the ICC result,
we judge that request worth making (item 3 of the programme), with the cyclic-versus-transitive
question as the specific target.
