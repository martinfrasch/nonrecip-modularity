# Note: plain-language verdict, and what it might say about Markov blankets

Saved 2026-08-25. Part 1 is the plain-language verdict on whether this work supports the
minAction/NWAP framework. Part 2 is speculative — hypothesis generation from an adjacent
system, not evidence about blankets.

---

# Part 1 — Does the framework stand supported?

**Short answer: partly — and the part that fails is the part that made it an *action* principle.**

## What was actually being claimed

The appeal of the idea was a rescue. Nonreciprocal forces have no potential energy function,
so you cannot write down an action and minimise it — that is the whole problem. The proposed
way out: split the interaction into two parts.

- a **symmetric** part ("we push each other equally"), which behaves like ordinary energy and
  *does* have an action
- an **antisymmetric** part ("I push you harder than you push me"), assumed to be *solenoidal*
  — stirring things in circles without changing what structures form

If that were true you keep your action principle. The symmetric part picks the structure by
energy minimisation; the antisymmetric part adds swirl on top without disturbing the answer.
The prediction followed: structure (modularity) frozen, circulation carrying the whole
nonreciprocal signature.

## What we found

**The antisymmetric part is what picks the structure.** Turning it up — with the symmetric part
held byte-for-byte identical — takes the system from 11 clusters to 137. Not a perturbation on
an energy landscape; the dominant term. The rescue does not work.

This had to be true in hindsight. Hara et al.'s whole result is that nonreciprocity *arrests
coarsening* — it changes what forms. A framework predicting that nonreciprocity cannot change
what forms was contradicting the experiment it was built to explain.

**The circulation is not there.** Looked for directly, twice, at two system sizes, against a
reference state where circulation is *provably* exactly zero. Nothing above noise.

**The proposed measurement failed backwards.** Newman modularity excess was supposed to be
elevated in the nonreciprocal case. It is *highest* in the reciprocal control and falls steadily
as nonreciprocity rises. Modularity is a broken ruler here — it reads ~0.86 whether the system
has 3 clusters or 33.

**The picture behind it does not hold either.** The energy-vs-connection-cost tradeoff needs the
system to settle on a preferred intermediate cluster size. Across a 4× range of system size it
does not — it phase-separates, the condensate growing in proportion to the box.

## What genuinely survived

**The decomposition itself was the right move.** Splitting the coupling into symmetric and
antisymmetric parts is exactly what this problem needed. It is what let us build a knob that
varies nonreciprocity and *nothing else* — which neither originally proposed knob did (the
field-strength parameter cancels out of the ratio algebraically; the size-ratio sweep, as the
paper parameterises it, varies packing instead).

With that knob, the central quantitative claim — **cluster statistics are set by the
antisymmetric-to-symmetric ratio** — is confirmed at p = 2×10⁻⁸, at two system sizes, with a
control the original experiment cannot perform. Real, publishable, and it came directly from the
framework's central structural idea.

**There is a better home for the antisymmetric part.** What it actually produces is
*dissipation* — measurable, positive, scaling as the square of the coupling. That quadratic form
is the signature of a different and well-established class of variational principle (Onsager's,
minimising energy change plus half the dissipation rate) rather than a least-action principle.

**Verdict: not an action principle for this system — but plausibly a dissipation principle, and
we measured the exact quantity such a principle would need.** A redirection, not a refutation.

## Caveats

One system, in simulation. It constrains what a general framework can claim about nonreciprocal
systems; it does not settle the framework's fate elsewhere. And the phase-separation result is
not finished — the two largest boxes are still slowly coarsening, continuations running. The
bias runs *against* the claim being tested, so a reversal is unlikely, but it is not closed.

## One-sentence version

The instinct to decompose the network into symmetric and antisymmetric sectors was right and is
validated; the story attached to the antisymmetric sector — that it circulates without
restructuring — is wrong, and what it does instead is dissipate, which points the framework
toward Onsager rather than toward least action.

---

# Part 2 — Clues about Markov-blanket systems?

**Speculative.** This system was never designed to test anything about blankets, and the mapping
from mechanical forces to conditional independence is not direct. What follows is hypothesis
generation, with measurements attached so the hypotheses are checkable.

## The useful thing here is a negative example

A Markov blanket requires *directed* influence: sensory states carry external→internal, active
states carry internal→external. A system whose couplings are perfectly reciprocal has a symmetric
influence graph and cannot support that asymmetry — every influence is matched by its reverse.
So **nonreciprocity is a necessary condition for a blanket.**

Our χ=0 state makes this concrete rather than merely definitional: detailed balance holds
*provably* (conservative central forces, uniform D/μ), entropy production is exactly zero, and
the stationary distribution is Boltzmann. There is no direction to anything. No blanket is
possible, not as a modelling choice but as a theorem.

Our χ>0 state clears that bar — directed pair influence, positive entropy production, broken
time-reversal symmetry. And yet **it has no blanket.** No internal/external partition, no
persistent boundary, no organised interface. It is a condensate shedding dimers.

That makes this system a clean example of **nonreciprocity without a blanket**, which is more
informative than another example of a system that has one. It isolates what directed coupling
buys you on its own — dissipation, restructuring, fragmentation — and what it conspicuously does
not. The gaps are candidate additional requirements.

## Gap 1 — the directed influence is incoherent

We measured the coherence of the antisymmetric forces inside a cluster:
C ≡ |Σ**p**| / Σ|**p**|. It satisfies C·√N ≈ constant from N=5 to N=80. The nonreciprocal
pushes add up **like a random walk at every scale** — they do not organise.

A blanket is precisely the opposite: the directed influences must be *coherently arranged*
around a closed surface, all the sensory channels pointing one way through it and all the active
channels the other. That is a highly ordered arrangement of exactly the quantity we measured to
be disordered.

**Candidate signature:** coherence of the antisymmetric coupling as a function of coarse-graining
scale. Unorganised nonreciprocity gives C ~ N^−1/2 (ours). A blanket-bearing system should show
C staying high up to the blanket's scale and dropping beyond it — the scale at which it drops
*is* the blanket's size. This is measurable in any model where the directed coupling is explicit.

## Gap 2 — no scale is selected

A blanket needs a characteristic size: the system must commit to how big the thing with an inside
and an outside is. Ours does not. Across a 4× range of system size it phase-separates — the
majority phase grows proportionally with the box (largest cluster ∝ N^0.75, lower-bounded by
incomplete coarsening, with the best-converged pair giving ∝ N^0.97) and the fragment population
is extensive (n_cl ∝ N^0.96, median size 3 at every box). There is no interior peak in the size
distribution at any scale tested.

The physical reason is plain: the model has attraction and nothing that opposes it at long range.
Pure attraction phase-separates. Systems that *do* select a scale have a competing long-range
term — screened repulsion, surface tension against internal pressure, a regulated division
mechanism, long-range inhibition à la Turing.

**Candidate requirement:** blanket formation needs a scale-selecting mechanism that is
independent of the nonreciprocity. Nonreciprocity supplies the direction; something else must
supply the size. Our system has the first and not the second, and gets no boundary.

## Gap 3 — irreversibility stays microscopic

This is the sharpest one, and it reframes our most frustrating null result.

We measured entropy production: positive, real, scaling as χ². The system is genuinely
irreversible. We *also* measured circulation in coarse observables — the signed area rate in
several observable planes, which is exactly zero under detailed balance — and found nothing
distinguishable from the equilibrium null, at two system sizes.

So: **irreversibility at the pair scale, none of it surviving coarse-graining.** Every microscopic
push is directed; none of that direction organises into a macroscopic cycle.

A Markov blanket is, thermodynamically, a *macroscopic* persistent cycle:
sensory → internal → active → external → sensory. If a blanket exists and is doing anything, the
directed flow must survive coarse-graining to the level of the partition. That is exactly the
signal we searched for and did not find.

**This suggests the circulation measure is a poor probe of nonreciprocity — our system is
strongly nonreciprocal and shows none — but may be a good probe of blanket-like organisation.**
It answers the question "does the directed structure survive coarse-graining?", which is closer
to what a blanket claim asserts than any structural measure is.

## A concrete suggestion for the dyad work

If the transfer-entropy analyses on maternal–fetal dyads are the empirical handle for blankets,
the parallel is direct. Entropy production and transfer entropy both measure directed,
irreversible structure — one thermodynamically, one informationally.

Our result in TE language would read: **strong pairwise directed influence, no directed influence
between coarse-grained partitions.** That is the fingerprint of nonreciprocity without a blanket.

So the discriminating measurement is not "is there directed influence" — there will be — but
**does directed influence persist between coarse partitions as you coarse-grain?** A blanket
predicts TE between the putative internal and external sets, conditioned on the blanket, falls to
zero (that is the definition), while TE *through* the blanket stays high and stays high under
coarse-graining. Unorganised directed coupling gives the opposite: micro-level TE everywhere,
washing out under coarse-graining.

## Where this could be wrong

- Markov blankets are a contested formal construct; whether physical systems literally have them
  versus have them assigned as a modelling convenience is actively disputed. Nothing here bears
  on that dispute.
- We measured coherence of *forces*. The mapping from mechanical coherence to conditional
  independence is an analogy, not a derivation.
- Our circulation null is low-powered. "No macroscopic circulation" here means "none detectable
  with these observables at this statistics", and the choice of observable planes was ours.
- The χ² dissipation scaling indicates *linear response* — the system responds passively and
  proportionally to being driven. Whether a blanket-bearing system should show departures from
  quadratic scaling (thresholds, saturation, indicating structured rather than passive response)
  is an appealing thought and entirely untested.

## Summary of the three candidate requirements

| | our colloids | a blanket would need | measurement |
|---|---|---|---|
| directed coupling | yes (χ>0) | yes | entropy production > 0 |
| coherent organisation of it | **no** (C ~ N^−1/2) | yes, up to the blanket scale | coherence vs scale |
| a selected size | **no** (phase separates) | yes | interior peak in size distribution |
| survives coarse-graining | **no** (no macro circulation) | yes | area rate / conditional TE between partitions |

Nonreciprocity gets you the first row and, on this evidence, none of the others by itself.
