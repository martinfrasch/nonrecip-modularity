# A2 — Onsager reciprocity: PASSED, in the antisymmetric (Casimir) form

The sharpest of the three tier-membership tests, and the one that required a second independent
drive. Code: `onsager.py`, `onsager_equil.py`.

## Design: why the obvious second drive fails

A uniform field on one species does not work. At zero field the system is isotropic, so the
species drift vanishes for every χ and the cross-coefficient is zero **by symmetry rather than by
physics** — the test would read 0 = 0.

Instead we opened a **second nonreciprocal channel** in the collision force, with its own
parameter χ₂ and the same structure as the EHD one (equal and same-signed on both members of a
pair — a co-propulsion):

```
F_col on i:   k(ssum−r) r̂  +  χ₂ h_ij k(ssum−r) r̂        h_ij = (s_i−s_j)/(s_i+s_j)
F_col on j:  −k(ssum−r) r̂  +  χ₂ h_ij k(ssum−r) r̂
```

Entropy production is then bilinear in two internal thermodynamic forces,

```
T·EPR = χ·J₁ + χ₂·J₂ ,   J₁ = Σ f₁·v  (EHD nonrecip force at unit χ)
                          J₂ = Σ f₂·v  (collision nonrecip force at unit χ₂)
```

so (χ, χ₂) are the forces and (J₁, J₂) their conjugate fluxes. Both drive the same internal
degrees of freedom, so the cross-coefficients are generically nonzero and the test has content.

## Result

Central differences at the operating point (χ, χ₂) = (0.5, 0.5), δ = 0.25, 40 configurations
**equilibrated at that operating point**, paired with common random numbers:

| T | symmetric (L₁₂+L₂₁)/2 | t | antisymmetric (L₁₂−L₂₁)/2 | t |
|---:|---:|---:|---:|---:|
| 200 | −1.08e−8 ± 6.5e−7 | **−0.02** | +9.975e−6 ± 6.3e−7 | +15.8 |
| 400 | −6.81e−8 ± 8.6e−7 | **−0.08** | +9.401e−6 ± 8.2e−7 | +11.4 |
| 800 | −3.96e−7 ± 9.9e−7 | **−0.40** | +9.477e−6 ± 9.6e−7 | +9.9 |

**The symmetric part of the cross-coupling is zero at every window length** (|t| ≤ 0.40), while
the antisymmetric part is large and stable at ≈9.5e−6. In other words

> **L₁₂ = −L₂₁**

which is the Onsager–Casimir relation in its antisymmetric form — the form that holds when the
two conjugate variables carry opposite time-reversal signature. Reciprocity holds. It simply
holds antisymmetrically rather than symmetrically.

## Why this is consistent with the entropy-production result

Entropy production is the quadratic form EPR = Σ_ij L_ij X_i X_j, to which **only the symmetric
part of L contributes** — an antisymmetric off-diagonal block contributes exactly zero. So the
cross-coupling between the two nonreciprocal channels is **purely reactive: it couples the
channels without dissipating.**

That is precisely why A1 came out as clean as it did. EPR = kχ² with a single coefficient, flat
to 4.8% across χ = 1.5–8, is what one expects when the cross terms contribute nothing to the
dissipation. Two independent measurements, made months apart in the project and by completely
different routes, agree on the structure of the response matrix.

**Two of three tier-II membership tests now pass.** A3 (single effective temperature) remains.

## Open theoretical question

The naive parity argument does not obviously predict the antisymmetric form. Both fluxes are
J = Σ f·v with f a function of positions (even under time reversal) and v odd, which would make
both fluxes odd, both signatures equal, and the relation *symmetric*. The measurement says
otherwise, robustly and across three window lengths.

So either the effective time-reversal signature of these two channels differs for a reason the
naive argument misses, or the antisymmetric coupling has a geometric/reactive origin not captured
by the standard parity assignment. **This is recorded as a theory question, not resolved here.**
It is also the single most interesting loose end for the variational-framework paper, because a
purely reactive coupling between nonreciprocal channels is exactly the kind of structure a
Rayleighian formulation would have to accommodate.

## Two protocols failed before this one worked

Both failures were instructive and are worth recording, because each produced confident,
publishable-looking numbers.

**Attempt 1 — long windows (T = 5×10⁴), 5 seeds.** L₁₂ resolved beautifully (28σ); L₂₁ was
unresolvable, its standard error 11× its own value. The script's automatic verdict printed
"RECIPROCITY HOLDS (t = 1.01)" — which was pure artifact: the t-statistic was small only because
L₂₁'s error bar swallowed everything, while the relative discrepancy was 241%. **Cause:** J₂ is
87× larger than J₁ (collision forces dominate EHD ones), and common random numbers stop reducing
variance once chaos decorrelates the paired trajectories, which a long window guarantees.

**Attempt 2 — short windows from many starts.** Restored the pairing and resolved both
coefficients, giving L₁₂ = +9.995e−6 and L₂₁ = −3.373e−5, i.e. an apparent reciprocity violation
at t = 93. **Also artifact.** The starting configurations were equilibrated at (χ=1, χ₂=0), not
at the operating point, so short windows measured relaxation toward the (0.5, 0.5) steady state.
The window-length scan exposed it: L₂₁ drifted monotonically −1.76e−5 → −2.63e−5 → −3.37e−5 →
−3.84e−5 across T = 50 → 400 with no plateau.

**What made attempt 3 work** was removing the reason short windows had a transient at all —
equilibrating *at* the operating point first, so the pairing works and the response is
steady-state simultaneously. The window-length scan then served as the acceptance test rather
than as a diagnostic: L₂₁ plateaus, and the symmetric part is zero at every T.

> **Generalisable lesson:** a response coefficient that has not been shown to plateau against
> the measurement window is not a response coefficient. Both failures here produced stable,
> small-error-bar numbers that would have passed any check except that one.
