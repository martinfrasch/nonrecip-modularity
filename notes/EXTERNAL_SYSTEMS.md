# Two candidate systems for extending the tests: Vicsek and Forkosh et al.

Assessment written 2026-09-02.

---

## 1. Vicsek model — useful, but not for the tier tests

**What it is.** Self-propelled particles at fixed speed v, each aligning its heading to the mean
heading of neighbours within radius r, plus angular noise η. Three control parameters: density,
noise amplitude, and v Δt / r. It has a well-characterised order–disorder transition with an
order parameter v_a = |Σ**v**_i| / (N v).

**Why it cannot support A1–A3.** Every tier-membership test in this work is calibrated against a
*provable equilibrium reference*. The entire method rests on χ = 0 giving conservative central
forces, uniform D/μ, exact detailed balance and identically zero current, so that any measured
current at χ > 0 is attributable to the antisymmetric coupling and not to the reference. **Vicsek
has no such limit.** Self-propulsion drives it out of equilibrium at every parameter value; there
is no setting of noise or alignment that recovers detailed balance. Without a calibrated zero
there is nothing to subtract, and the entropy-production estimator — which already required a
paired same-configuration control to be usable at all — has no anchor.

Note also that standard Vicsek alignment is **reciprocal**: each particle averages over its
metric neighbourhood, mutually. Nonreciprocity has to be added, most naturally as a *vision cone*
(a particle aligns only to neighbours ahead of it), which makes the interaction non-mutual because
i may see j while j does not see i.

**What it is genuinely good for, and it is important.** Our colloid-versus-cilia comparison
(Section 3.11) confounds two distinct contrasts: biological versus synthetic, and *collectively
ordered* versus *spatially distributed and uncoordinated*. A referee raised exactly this. Vicsek
resolves the confound because it is **synthetic and collectively ordered**.

The sharpest version uses **chiral Vicsek** — self-propelled particles with an intrinsic turning
rate, i.e. circle swimmers — which produces a collective *limit cycle*, structurally the same
object as the ciliary beat. The test is then:

> Apply the identical signed-area-rate estimator, at matched level of description, to a chiral
> Vicsek flock.

- If it shows circulation like cilia, the discriminating property is **limit-cycle collective
  organisation**, not biology, and our Section 3.11 claim must be reframed accordingly.
- If it shows none, the biological reading survives a serious attempt to break it.

This is a control that can *weaken* our own headline, which is precisely why it should be run. It
is also cheap: Vicsek is a few lines and runs in minutes at N = 10⁴.

**Verdict:** adopt Vicsek as a coarse-graining control, not as a tier-test system. Add a chiral
variant as the synthetic limit-cycle positive control that Section 3.11 currently lacks.

---

## 2. Forkosh et al., *Nat. Neurosci.* **22**, 2023 (2019) — the mathematics does not connect; the system does

**The framework.** "Identity domains" are found by seeking directions in behaviour space with
maximum consistency and discriminative power — that is, by maximising the ratio of
between-individual to within-individual variance, a generalised eigenvalue problem of
Fisher-discriminant type. It is a variance decomposition, not a variational principle in the
physical sense: there is no action, no path measure, no dissipation, no time-reversal operation,
and no conserved or produced quantity. The formal resemblance to NWAP — both extremise something —
is superficial, and we would not claim a connection at the level of the mathematics. Doing so
would be exactly the kind of surface-level analogy that this project has otherwise been careful
to avoid.

**Where the real connection lies.** Group-housed mice are a *nonreciprocal interacting many-body
system*, and dominance relations are the canonical biological instance of broken action–reaction
symmetry: if A chases B, B does not chase A. The published behaviour table records this directly.
Of 67 readouts, 31 are social, and several come in explicitly **directed pairs**:

| forward | reverse |
|---|---|
| `AggressiveChaseRate` | `AggressiveEscapeRate` |
| `FollowRate` | `BeingFollowedRate` |
| `FractionOfChasesPerContact` | `FractionOfEscapesPerContact` |
| `DiffBetweenApproachesAndChases` | — (already an asymmetry) |

That is the antisymmetric sector of a social interaction, measured.

**The limitation, and it is decisive for a full test.** The public table
([github.com/OrenForkosh/IdentityDomains](https://github.com/OrenForkosh/IdentityDomains), MIT
licence) is 848 rows × 67 columns — **aggregated per mouse per day**, four days, with group and
mouse identifiers. It gives each animal's rate of chasing and of being chased, but **not who
chased whom**. The directed adjacency G_ij is therefore not recoverable from it, and neither a
directed-network analysis nor a transfer-entropy measurement between named individuals is
possible without the raw tracking data.

**What is possible with the public data:**

- A per-individual asymmetry index (chase minus escape, follow minus being-followed), and its
  consistency across days — testable, since each animal appears on four days.
- Whether that asymmetry index aligns with the identity domains, i.e. whether "dominance" emerges
  as a principal axis of the variance decomposition. This is a real question and cheap to answer.
- Group-level asymmetry statistics: whether the distribution of individual asymmetries is
  consistent with a transitive hierarchy or with cyclic (rock–paper–scissors) structure. Cyclic
  dominance would be a genuine social analogue of circulation.

**What would need the raw data:** the directed interaction matrix, from which one could measure
whether irreversibility survives coarse-graining in a social system — the same test applied to
colloids and cilia, at a third and very different scale. That would be the substantive extension,
and it requires contacting the authors.

**Verdict:** do not claim a mathematical connection. Do treat the system as a candidate third
test-bed, note that the public data supports only the individual-level asymmetry analysis, and
approach the authors if the directed-network test is worth pursuing.

---

## Recommended order

1. **Chiral Vicsek circulation control.** Cheap, fast, and it directly addresses a referee
   objection to Section 3.11 that no amount of additional colloid work can answer. It is also the
   experiment most likely to change our conclusion.
2. **Forkosh asymmetry index against the identity domains.** A few hours' analysis on already-open
   data; answers whether directed social asymmetry is a principal axis of behavioural individuality.
3. **Raw-tracking request**, only if (2) is suggestive.

Neither system can substitute for the colloid model in the tier tests, because neither has an
equilibrium reference. That remains the distinguishing methodological asset of the present work.
