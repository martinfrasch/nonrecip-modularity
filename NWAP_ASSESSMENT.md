# The Hara et al. system as a NWAP test: what 82 runs say

Assessment of the three claims in the original review of Hara et al. (PRL 137, 068302, 2026)
against the simulation evidence. Sources: `FINAL_RESULTS.md`, `AUDIT.md`, `COARSENING.md`.

**Headline: claim 2's quantitative core is confirmed and is the strongest result here.
Claim 1's proposed empirical test fails, in the wrong direction. Claim 3 gains a sharper
and measurable form. But the mechanism attributed to the antisymmetric sector — solenoidal,
circulating, structure-preserving — is wrong on all three counts, and that changes what kind
of variational object can accommodate this system.**

---

## Claim 1 — "arrested coarsening is your modularity signature in a dish"

**Proposed test:** compute Newman modularity of the contact network against an equilibrium
(monodisperse) null; NWAP predicts *sustained modularity excess in the bidisperse case*.

**Result: fails, and the sign is inverted.** Modularity excess over the degree-preserving
null, paper scale:

| | mono (equilibrium null) | χ=0 | χ=0.25 | χ=0.5 | χ=0.75 | χ=1 | χ=1.5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ΔQ vs null | **0.462** | 0.440 | 0.423 | 0.409 | 0.395 | 0.389 | **0.383** |

Modularity excess falls monotonically as nonreciprocity rises, and is **highest** in the
equilibrium null. The prediction is not merely unsupported; the ordering is reversed.

**Why:** Q is degenerate on these graphs. Across all conditions run, cluster count spans 40×
while Q spans 14%, and the two are uncorrelated (ρ=+0.31, p=0.28). A 3-cluster and a
33-cluster state are indistinguishable in Q. Louvain's resolution limit merges anything below
~√(2E) nodes, and this survives a 4× increase in system size. Newman modularity cannot carry
this signature.

**The deeper problem — there is no intermediate characteristic scale.** The E-vs-C tradeoff
picture requires the system to select a finite cluster size: an interior peak in the
size distribution. Late-window distributions (paper scale):

| χ | clusters/snapshot | median size | fraction of particles in system-scale clusters | log-size histogram |
|---|---:|---:|---:|---|
| 0 | 10.5 | 214 | 0.72 | 2,0,0,16,43,43,1,11 ← **interior peak** |
| 0.25 | 3.9 | 507 | 1.00 | condensed |
| 1 | 98.5 | **3** | 0.84 | 782,200,72,13,5,0,0,11 ← monotonic |
| 1.5 | 135.5 | **3** | 0.82 | 1186,204,68,16,5,0,0,11 ← monotonic |

The nonreciprocal steady state is **one system-scale condensate holding ~82% of particles,
continuously shedding a monotonically-decaying spray of dimers and trimers** (median cluster
size 3). That is not a population of finite reconfiguring clusters at a selected scale.

Ironically the *only* condition with an interior peak is χ=0 — and that peak is set by
diffusion-limited aggregation kinetics (`COARSENING.md`), not by any E-vs-C tradeoff.

**Caveat, and it matters:** the condensate is system-scale in a box of 24. If the true
characteristic scale exceeds the box, a monotonic distribution is what a too-small system
would show regardless. Testing claim 1 properly requires L ≫ the condensate size; that is
untested and is the single most valuable follow-up run.

**Verdict:** the proposed reanalysis of Hara's published data via modularity excess should not
be attempted — we ran it and it fails. Cluster-size statistics, not modularity, carry the
signature.

---

## Claim 2 — "nonreciprocity is a directed graph; decompose G_ij into symmetric and antisymmetric parts"

**This is right, and it is the core contribution.** The decomposition is exactly what the
system needed, and implementing it is what made every subsequent result possible.

Operationally: χ scales the antisymmetric sector while the symmetric part stays *bit-identical*
(verified — an equal-radius pair follows identical trajectories at χ=0 and χ=1 to 1e-14).
χ=0 is exactly reciprocal, χ=1 is the original Hara force.

### The central quantitative prediction is confirmed

*"The cluster-size distribution and reorganization rate should be set by the ratio of
antisymmetric to symmetric coupling strength."*

| quantity vs χ | Spearman ρ | p |
|---|---:|---:|
| cluster count | **+0.931** | 2.1e-08 |
| σ²_v (activity) | **+0.975** | 7.1e-12 |

n_cl goes 11.1 → 137.3 (12.4×, p=9.4e-10) with everything else held fixed, at two system
sizes. **Cluster statistics are set by the antisymmetric/symmetric ratio.** This is the one
NWAP prediction that survives contact with the data, and it survives cleanly.

"Reorganization rate" is weaker: σ²_v tracks χ, but edge turnover — the natural network
reorganization observable — is **non-monotonic** at paper scale, peaking at χ=0.75. It is also
morphology-correlated (ρ=+0.795 with n_cl) and ~74% of it survives at χ=0. Do not use it.

### But the attributed mechanism is wrong in three ways

1. **"The antisymmetric part cannot alter a static structural observable."** False. n_cl is a
   single-snapshot observable and it moves 12.4×. The premise — antisymmetric ⇒ solenoidal ⇒
   density-preserving — holds only when ∇·(ρ_eq **v**) = 0, a special condition generic
   nonreciprocal couplings do not satisfy. It also contradicts Hara et al.'s own result:
   arrested coarsening *is* a nonreciprocity-induced change in steady-state structure.

2. **"The antisymmetric part generates the circulating component."** Not detected. The signed
   area rate in observable planes is exactly zero under detailed balance; measured in four
   planes at every χ at both scales, everything is indistinguishable from the χ=0 equilibrium
   null (all p ≥ 0.38). "Circulating partition" has no empirical support here.

3. **What the antisymmetric sector actually produces is dissipation.** Entropy production —
   which unlike modularity, turnover and σ²_v is *exactly zero* under detailed balance — is
   positive for χ ≳ 0.75 and, at fixed configuration, **scales as χ²** (net/χ² flat to 1.27×
   and 1.61× on two of three seeds). The dynamics are genuinely irreversible. The current is
   real; it simply does not appear as circulation in any network observable.

### Correction to the proposed experimental knob

"Experimentally tunable via the size ratio and field frequency" — with one caveat that
matters for anyone repeating the simulation:

- **α̂ (field-strength analogue) cancels exactly** from the antisymmetric/symmetric ratio.
  Both sectors scale linearly in it, so it tunes overall EHD strength and nothing else
  (verified identical to 6 d.p. across α̂ ∈ {0.003…0.015}). Sweeping it produced a trend
  *opposite* to prediction, for reasons unrelated to reciprocity.
- **The size ratio is a valid experimental knob but not a valid simulation knob** as the paper
  parameterises it. In experiment a particle's EHD radius and steric radius are the same
  physical size, so changing the size ratio changes both and does tune the nonreciprocity.
  In the paper's S5/S8 protocol the EHD radii are *pinned* while steric radii vary, which
  varies packing only. Our s_II/s_I sweep therefore tested nothing about reciprocity.

The knob that works in simulation is χ, or equivalently the EHD radius contrast l_S/l_L.

---

## Claim 3 — the biological bridge

The results sharpen this into something measurable. **There are two physically distinct kinds
of arrest in this system, and they are indistinguishable structurally but separated cleanly by
dissipation:**

| | χ=0 | χ ≳ 0.75 |
|---|---|---|
| state | arrested | arrested |
| mechanism | kinetic — 5–7 kT bonds, cluster diffusion dying as 1/N | active — continuous fragmentation/reformation |
| entropy production | **exactly zero** (detailed balance holds provably) | **positive, ∝ χ²** |
| stability | a dead gel; the condensed state is equally stable and it simply cannot reach it | maintained only while driven |

This is precisely the frozen-aggregate versus homeostasis contrast, and it is now
operational: the discriminator is entropy production, not structure. A frozen aggregate and a
dynamically maintained one can look similar in any static network measure; they differ
absolutely in dissipation. For the maternal–fetal and directed-synapse analogies, that
suggests the transferable quantity is an irreversibility measure, not a modularity measure.

One caution for the analogy: the χ=0 gel is *also* "arrested", so "arrested coarsening" alone
does not indicate an active state. The paper's monodisperse reference behaves the same way
(lcf 0.30, n_cl 11.2, statistically identical to χ=0). Any claim of the form "arrest indicates
sustained organisation" needs the dissipation measurement to back it.

---

## Verdict on the framework

The original verdict — *"a framework whose weights live on an undirected network fails here;
NWAP survives iff you formalise the directed-G extension"* — is right about the necessity of
the extension and wrong about what the extension buys.

**The directed-G decomposition is necessary and it works.** Its quantitative prediction —
cluster statistics set by the antisymmetric/symmetric ratio — is confirmed at p=2e-8 across
two system sizes. That is a real, falsifiable, reproduced result and it is publishable.

**But the antisymmetric sector cannot be assigned the role the framework gives it.** You
cannot partition "symmetric → energy/structure, antisymmetric → circulation/information flux".
Both sectors set the structure, and the antisymmetric one dominates it. What the antisymmetric
sector uniquely contributes is a **positive quadratic dissipation**: at fixed configuration
EPR = k·χ².

That functional form is the useful finding. A positive quadratic form in the antisymmetric
coupling is Onsager/Rayleigh structure, not action-extremum structure. So the constructive
recommendation:

- Put the antisymmetric sector in a **dissipation functional** (a Rayleighian, or a
  steady-state entropy-production principle), not in an action whose extremum is supposed to
  leave the stationary density invariant. The measured EPR ∝ χ² is exactly the object such a
  formulation predicts.
- The **I (information flux) term needs a different observable.** Both candidates tried here —
  modularity excess and partition/edge turnover — failed: the first has the wrong sign, the
  second is non-monotonic and largely a morphology readout. Entropy production is the
  measurable that behaves.
- The **E-vs-C tradeoff needs an interior peak in the cluster-size distribution to exist.**
  At accessible system size it does not. Establish it at L ≫ condensate size before building
  on it.

## For the proposed follow-up paper

The empirical hook as proposed — modularity-excess reanalysis of Hara's published data —
should be dropped; it was tested and fails with inverted sign.

The hook that works: **the χ decomposition as a controlled knob, with n_cl(χ) and EPR(χ).**
It is a clean directed-graph result, it isolates the antisymmetric sector with a
same-particle reciprocal control that Hara et al. do not have, and it yields both a confirmed
structural prediction and a quadratic dissipation law. The two-kinds-of-arrest distinction,
with entropy production as the discriminator, is the natural bridge to the biological claim.
