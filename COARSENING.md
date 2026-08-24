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

## What weak activity does: size-independent cluster propulsion

Thermal cluster mobility dies as 1/N. Nonreciprocal propulsion does not, because a cluster's
net driving force grows with its number of nonreciprocal pairs while its drag grows with N.
Measured cluster-COM speed by size:

| cluster size | χ=0 (thermal) | χ=0.75 | χ=1 |
|---|---:|---:|---:|
| 4–12 | — | 2.67e-4 | 2.86e-4 |
| 12–40 | — | 2.59e-4 | 2.95e-4 |
| 40–120 | **2.44e-5** | **2.64e-4** | **2.62e-4** |
| 120–400 | 1.24e-5 | 6.14e-5 | — |
| 400+ | 7.42e-6 | 7.89e-6 | 6.79e-6 |

Two things to read off:

- Active speed is **flat from N=4 to N=120** — a 30× size range at constant speed, exactly the
  size-independent propulsion the argument predicts. At matched size (40–120) active clusters
  move **10.7× faster** than thermal ones. That is the unjamming.
- Above ~120 particles it crosses back to thermal-like scaling and the largest clusters move
  identically regardless of χ. Propulsion directions decorrelate inside a big cluster, so net
  force grows as √N rather than N and velocity falls as N^−0.5 again. This is also why the
  condensate, once formed, stays put.

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

- The crossover size (~120 particles) where propulsion decorrelates is measured, not derived.
  It should follow from the persistence length of propulsion-direction correlations within a
  cluster, which is not computed here.
- χ=0.5 is the only condition still evolving at t=4e5 (+18% in n_cl over the last quarter),
  sitting right at the crossover. Its steady state is not established.
