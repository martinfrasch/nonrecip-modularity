# Experiment spec: does nonreciprocity move currents without moving structure?

Status: **RUN 2026-08-21.** Results and scoring against the pre-committed thresholds below
are in `RESULTS.md` §5. Outcome: H1 pass, H3 pass, H2 fail-as-written (decoupling confirmed
at 16x separation, but the absolute turnover threshold was set too high), H4 partial --
reciprocity accounts for ~half the mono/bi turnover gap and essentially all of the sigma^2_v gap.

Written after the alpha and s_II/s_I sweeps were found to be mis-specified instruments
(see `RESULTS.md` section 4).

---

## 1. Why the previous sweeps could not test this

The claim under test is that the antisymmetric part of the coupling is solenoidal: it moves
probability currents while leaving static structure alone. Testing that requires a knob on
the **antisymmetric/symmetric ratio**. Neither prior sweep was one.

- **α̂** multiplies the symmetric *and* antisymmetric parts equally, so it cancels exactly
  from their ratio (verified: ratio identical to 6 d.p. across α̂ ∈ {0.003…0.015}). α̂ is a
  coupling-**strength** knob. Raising it condensed the system (n_cl 40 → 9), which suppressed
  edge turnover — giving a trend *opposite* to the prediction for reasons unrelated to it.
- **s_II/s_I** moves steric radii with EHD radii pinned, so it varies packing, not reciprocity.

Both knobs move structure and currents together, which is exactly what was observed.

## 2. The knob

Split the EHD force into its reciprocal and nonreciprocal parts and scale only the latter.
With g_k ≡ α̂·l_k⁴/(r²+l_k²)^(5/2) and ḡ ≡ (g_i+g_j)/2:

```
c_i = (1−χ)·ḡ + χ·g_i          F_i = −c_i·r̂
c_j = (1−χ)·ḡ + χ·g_j          F_j = +c_j·r̂

symmetric part   (c_i+c_j)/2 = ḡ            ← independent of χ, identically
antisymmetric    (c_j−c_i)/2 = χ·(g_j−g_i)/2 ← linear in χ
```

- **χ = 0** — Newton's third law holds exactly. A reciprocal control with *the same particles*.
- **χ = 1** — recovers the original Hara et al. force exactly.
- **χ = 1.5** — extrapolation past the physical model; exploratory, flagged as such.

This is the only construction here that varies reciprocity while holding the energetic
backbone bit-for-bit fixed, for every pair type, every separation, every size.

Implemented in `simulate.py:_run` (5 lines) and validated end-to-end by
`tests/test_reciprocity.py` against the real kernel:

| check | result |
|---|---|
| χ=0 net pair force (isolated noise-free pair, COM drift) | **0.000e+00** — exactly reciprocal |
| χ=1 COM drift | 3.95e-2 — self-propelling |
| χ=0.5 drift vs χ=1 | 1.975e-2 = exactly ½ — linear in χ |
| equal-radius pair, χ=0 vs χ=1 trajectories | identical to 1e-14 — symmetric part invariant |

## 3. It also fixes a confound in the baseline

The published baseline compares **monodisperse vs bidisperse**. Those differ in *two* things
at once: reciprocity (equal particles are automatically reciprocal) and size polydispersity
(different packing, different cluster morphology, different N at matched area fraction).
The 2.4× turnover difference cannot be attributed to reciprocity alone.

χ=0 bidisperse vs χ=1 bidisperse is the missing control: **identical particles, identical
packing, identical symmetric coupling, differing only in reciprocity.** Whatever separates
those two arms is attributable to nonreciprocity and nothing else.

## 4. Design

| | |
|---|---|
| Factor | χ ∈ {0, 0.25, 0.5, 0.75, 1.0, 1.5} |
| Seeds | 5 per level (independent RNG streams; seeds 1–5) |
| Runs | 30 |
| Held fixed | α̂=0.005, s_II/s_I=0.667, l_L=1/6, l_S=1/9, N=1000, box 12, t̂=1e5, dt̂=0.05, 100 snapshots |
| Analysis window | late half (`--burn 0.5`) |
| Primary outcome | edge Jaccard turnover (algorithm-free) |
| Secondary | Q, ΔQ vs degree-preserving null, σ²_v, n_cl, lcf, ARI vs its own noise floor |

Command:
```bash
python simulate.py --sweep chi --workers 8     # ~55 min on 8 P-cores
python analyze.py
python figures.py --chi                        # -> chi_test.png
```

**Blocking on seeds:** seed k uses the same initial configuration across all six χ levels
(`np.random.default_rng(seed)` is drawn before the dynamics), so χ is compared within
matched initial conditions. Analyse as a paired design where possible.

## 5. Predictions, stated so they can fail

**H1 (primary).** Edge turnover increases monotonically in χ.
Test: Spearman ρ between χ and per-run edge turnover over all 30 runs.
Pass: ρ > 0 with p < 0.01.

**H2 (the actual claim — decoupling).** Structure moves far less than currents.
Test: relative change from χ=0 to χ=1, R_X ≡ |X(1) − X(0)| / X(0).
Pass: R_Q < 0.1 **and** R_turnover > 0.5, i.e. at least a 5× separation.
This is the "frozen Q, circulating partition" claim in falsifiable form.

**H3 (control).** σ²_v increases with χ, confirming the antisymmetric part injects current.

**H4 (null recovery).** At χ=0 the bidisperse system's turnover should fall toward the
monodisperse baseline value (0.164). If χ=0 still shows ~0.40 turnover, then the mono/bi
difference was driven by polydispersity, not reciprocity — which would **falsify the
project's headline interpretation** while leaving the measurement intact. This is the
single most informative arm and the reason χ=0 must be run.

**Pre-committed failure readings:**

| observation | reading |
|---|---|
| turnover flat in χ | nonreciprocity does not drive edge currents; headline claim dies |
| Q moves as much as turnover | no decoupling; "frozen Q" is wrong |
| χ=0 turnover ≈ χ=1 turnover | mono/bi baseline gap was polydispersity, not reciprocity |
| turnover rises, Q flat | **prediction confirmed** — first clean support in this project |

## 6. Power

Baseline within-condition seed scatter for edge turnover: sd ≈ 0.019, SEM ≈ 0.008 at n=5.
The χ=0→1 contrast is expected to span roughly 0.16 → 0.40 if the effect is real — ~12 SEM.
Five seeds is ample for the primary contrast; the 30-run Spearman is the conservative test.
Q's SEM is ~0.002, so H2's 10% threshold on Q (0.086 in absolute terms) is far above noise —
meaning H2 fails loudly if Q really does move.

## 7. Known limits

- χ>1 is outside the derived electrohydrodynamic model; treat χ=1.5 as trend evidence only.
- t̂=1e5 at box 12 is the pilot scale, not the paper's steady-state scale (t̂≥4e5, box 24).
  Absolute cluster sizes will not match the paper; the χ *contrast* is the object here, and
  all arms share the scale. A confirmed result should be repeated at paper scale before
  publication.
- Pairwise forces only — no many-body hydrodynamics (inherited from the paper, SI S6).
- Louvain's noise floor is condition-dependent; that is why turnover, not ARI, is primary.

## 8. If it confirms

The natural follow-up is the EHD radius-contrast sweep (l_S from 1/6 → 1/12 at fixed steric
radii), which varies the ratio through the *physical* channel rather than a synthetic mixing
parameter. Note it confounds S–S coupling strength (g ∝ l⁴), so it needs per-condition α̂
compensation to hold the symmetric part fixed — messier than χ, which is why χ is the
primary instrument and the radius sweep the corroboration.
