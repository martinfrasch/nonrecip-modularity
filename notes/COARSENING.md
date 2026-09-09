# Why cluster count dips at χ=0.25

The cluster-count non-monotonicity — n_cl falls from 11.1 at χ=0 to **2.2** at χ=0.25 before
climbing to 137 at χ=1.5 — reproduces at both system sizes and was flagged as unexplained in
`FINAL_RESULTS.md`. It has a mechanism, and it changes how the reciprocal reference state
should be read.

**It is not a dip in coarsening. It is coarsening finally completing.**

## The reciprocal state is a kinetically arrested gel, not an equilibrium reference

At χ=0 the system plateaus at largest-cluster-fraction **0.30** and stays there from
t≈160k onward — ten mid-sized clusters, no dominant one. That is not phase separation, and
the system is not still on its way there:

| | χ=0 | χ=0.25 |
|---|---:|---:|
| lcf at t=400k | 0.301 | **0.910** |
| n_cl at t=400k | 10.3 | **2.0** |
| slope over last quarter | −0.7% | −4.5% (both steady) |

Three independent lines of evidence say χ=0 is arrested rather than equilibrated:

**1. Bonds are effectively irreversible.** Integrating the χ=0 symmetric pair force gives well
depths of **−7.0 kT (L-L), −6.8 kT (L-S), −4.8 kT (S-S)** against kT_eff = σ²/2 = 3.38e-6.
Contacts essentially never break thermally, so aggregation is diffusion-limited and the
resulting network cannot rearrange.

**2. Large clusters cannot find each other.** A rigid cluster of n particles has D ≈ D₁/n:

| n | 1 | 10 | 100 | 1000 |
|---|---:|---:|---:|---:|
| rms displacement over a full t=4e5 run | 5.70 | 1.80 | 0.57 | **0.18** |

In a box of 24, a 1000-particle cluster moves less than one particle diameter in an entire
run. Coalescence stalls for purely kinetic reasons.

**3. The condensed state is stable at χ=0 — decisive.** Take the condensed χ=0.25
configuration and evolve it at χ=0 with no activity at all:

| transition | lcf start → end | n_cl start → end |
|---|---|---|
| **condensed → χ=0** | 0.910 → **0.910** | 2.0 → **2.0** |
| condensed → χ=0.25 (ctrl) | 0.910 → 0.910 | 2.0 → 2.0 |
| dispersed → χ=0 (ctrl) | 0.301 → 0.301 | 10.3 → 10.0 |
| **dispersed → χ=0.25** | 0.301 → **0.321** | 10.3 → **8.7** |

Per-seed trajectories for condensed→χ=0 are dead flat: `0.763→0.763`, `1.000→1.000`,
`0.966→0.967`. The condensate is perfectly stable without any activity, so it is
equilibrium-accessible — χ=0 simply cannot reach it from a dispersed start. Meanwhile the
same dispersed configuration is frozen at χ=0 but coarsens at χ=0.25.

That is kinetic arrest, not a distinct nonequilibrium steady state.

## What weak activity does: a uniform ~10× mobility boost (corrected)

**This section corrects an earlier claim.** I originally reported that active cluster speed is
*size-independent* up to ~120 particles, and proposed deriving that crossover from a
propulsion-correlation length. Measuring the driving force directly shows there is no such
crossover, and the size-independence was a measurement artifact.

### The propulsions add incoherently at every size

A cluster's drive is the vector sum of its pair propulsions **p**_ij = −χ(G_i−G_j)/2·**r**_ij
(identical on both members of a pair). Coherence C ≡ |Σ**p**| / Σ|**p**| distinguishes the
two regimes: C constant ⇒ aligned drive, |F|~N, speed independent of size; C ~ N^−1/2 ⇒
random-walk drive, |F|~√N, speed ~ N^−1/2.

| cluster size | ~5 | ~13 | ~34 | ~80 | ~3316 |
|---|---:|---:|---:|---:|---:|
| C | 0.801 | 0.488 | 0.339 | 0.287 | 0.015 |
| **C·√N** | **1.79** | **1.76** | **1.97** | **2.57** | 0.88 |

C·√N is flat across the whole small-to-medium range: **the pair propulsions add like a random
walk at every size.** There is no coherently-driven domain and therefore no coherence
crossover to derive.

### Consequently the drift scales exactly like thermal motion

Instantaneous drift speed from the force, v = |F_net| / Σs_i:

| N | ~5 | ~13 | ~34 | ~80 | ~3316 |
|---|---:|---:|---:|---:|---:|
| v_inst (from force) | 5.68e-4 | 4.14e-4 | 2.75e-4 | 2.33e-4 | 1.11e-5 |
| v_meas (Δt=2000 displacement) | 3.09e-4 | 2.99e-4 | 2.41e-4 | 2.37e-4 | 1.30e-5 |

Overall fit: **v_inst ~ N^−0.50** — the same exponent as thermal cluster diffusion.

### Where the artifact came from

v_inst and v_meas agree for large clusters (1.11e-5 vs 1.30e-5) but diverge for small ones
(5.68e-4 vs 3.09e-4). Over a Δt=2000 snapshot interval a small cluster reorients many times,
so its *net displacement* underestimates its instantaneous speed; a large cluster moves
ballistically over the same window and does not. That size-dependent suppression is what
flattened the measured curve and produced an apparent plateau. Displacement over a coarse
interval is not a speed — the same class of error `AUDIT.md` documents for σ²_v.

### What survives

Activity changes the **amplitude**, not the scaling. At matched cluster size (40–120
particles) active clusters move **10.7× faster** than thermal ones, and that ratio holds
across sizes because both scale as N^−1/2. A ~10× mobility boost is more than enough to
restart stalled coalescence, so the unjamming account below is unaffected — only its
mechanism is corrected from "size-independent propulsion" to "uniform mobility enhancement".

## The full picture

Activity has two competing effects, both increasing with χ:

1. **mobilises clusters** → coalescence resumes → *fewer, larger* clusters
2. **tears clusters apart** → fission → *more, smaller* clusters

At χ=0.25 mobility wins and the system completes the phase separation that χ=0 could not
reach. By χ≳0.75 fission dominates and cluster count climbs monotonically to 137. The dip is
the crossover between the two.

## Consequence for the project's headline result

**The reciprocal reference state is not an equilibrium control — it is a gel.** Comparisons
against χ=0 (and the paper's monodisperse comparison, which behaves the same way: lcf 0.30,
n_cl 11.2) partly measure *unjamming* rather than nonreciprocity as such.

This does not overturn the arrested-coarsening result — n_cl still rises 12.4× from χ=0 to
χ=1.5, and from χ=0.25 (the true coarsened state, n_cl=2) it rises **~60×**, which is a
larger effect measured against the right baseline. But it reframes it: relative to the true
phase-separated state, **both** ends are arrested — χ=0 kinetically, high χ actively. What
nonreciprocity sets is the arrest *mechanism* and the arrest *scale*, not arrest versus
coarsening.

Anyone quoting "nonreciprocity arrests coarsening" from these runs should state which
baseline they mean.

## Reproduce

```bash
python arrest_test.py --workers 8          # config-swap test -> arrest_test.csv
python figures.py --coarsening             # -> coarsening.png
```

## Open

- The ~120-particle crossover does not exist (see the corrected section above); it was an
  artifact of coarse-grained displacement. The real scaling is N^−1/2 at all sizes.
- The local exponent does steepen for the very largest clusters (−0.82 in the top bin versus
  −0.50 overall). Those clusters hold 83% of all particles and span 23 of 24 box slabs, so
  finite-size pinning is the likely cause, but this is not established.
- χ=0.5 was flagged as "still evolving at t=4e5 (+18% over the last quarter)". **Resolved:
  it is steady.** Continuing all three seeds to t=8e5 gives n_cl = 17.3 ± 0.5 in the late half
  of the first window and 18.8 ± 0.9 in the late half of the second (Welch p=0.23). The
  apparent trend was a quarter-window slope fit against ~19% intra-run fluctuation.

  Fluctuations are strongly non-uniform across the sweep and peak in the crossover region:

  | χ | 0 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 |
  |---|---:|---:|---:|---:|---:|---:|
  | sd/mean of n_cl (late window) | 7.2% | 29.7% | 18.6% | 11.9% | 9.3% | 7.1% |

  Minimal deep in either regime, maximal where activity-driven coalescence and
  activity-driven fission balance. (χ=0.25's figure is inflated by its small mean — n_cl≈2,
  so ±1 cluster is 50% — but the qualitative pattern holds in absolute terms too.) Any
  convergence test in the crossover region needs a window long enough to average over this;
  a quarter-window slope will read as a trend when there is none.
