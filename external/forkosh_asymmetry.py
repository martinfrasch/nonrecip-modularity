"""Is directed social asymmetry a stable individual trait, and does it behave like a
nonreciprocal coupling?

Group-housed mice are a nonreciprocal interacting system: if A chases B, B does not chase A.
The published behaviour table (Forkosh et al., Nat. Neurosci. 22, 2023, 2019; data at
github.com/OrenForkosh/IdentityDomains, MIT) records both directions of several social
interactions per animal per day, which is the antisymmetric sector of a social coupling,
measured.

What the public table does NOT contain is who interacted with whom -- it is aggregated per mouse
per day -- so the directed adjacency G_ij is not recoverable and no transfer-entropy or
directed-network measurement is possible from it. What is possible:

  1. a conservation check that validates the interpretation (within a closed group, every chase
     is someone's escape, so the asymmetries must sum to zero)
  2. whether the per-individual asymmetry is a stable trait across days, or day-to-day noise
  3. how steep the resulting hierarchy is, and whether steepness is itself consistent
"""
from __future__ import annotations
import numpy as np, pandas as pd

PAIRS = [("AggressiveChaseRate", "AggressiveEscapeRate", "aggression"),
         ("FollowRate", "BeingFollowedRate", "following"),
         ("FractionOfChasesPerContact", "FractionOfEscapesPerContact", "chase/contact"),
         ("FractionOfFollowsPerContact", "FractionOfBeingFollowedPerContact", "follow/contact")]


def load():
    d = pd.read_csv("external/data_forkosh/behaviors_table.csv")
    for f, r, name in PAIRS:
        d[f"asym_{name}"] = d[f] - d[r]
    return d


def icc(df, value, unit="MouseNumber"):
    """Intraclass correlation: between-individual variance as a fraction of the total.
    This is the same between-over-within variance logic the identity-domains framework
    maximises, applied to one pre-specified axis instead of searching for the best one."""
    g = df.groupby(unit)[value]
    k = g.size().mean()
    msb = g.mean().var(ddof=1) * k
    msw = g.var(ddof=1).mean()
    return max(0.0, (msb - msw) / (msb + (k - 1) * msw))


if __name__ == "__main__":
    d = load()
    print(f"=== {d.MouseNumber.nunique()} mice, {d.GroupNumber.nunique()} groups, "
          f"{d.Day.nunique()} days, {len(d)} rows ===\n")

    print("1. Conservation check: within a closed group every chase is another animal's escape,")
    print("   so the per-group sum of asymmetries must vanish.\n")
    for f, r, name in PAIRS:
        s = d.groupby(["GroupNumber", "Day"])[f"asym_{name}"].sum()
        sc = d[f"asym_{name}"].abs().mean()
        print(f"   {name:14s} mean |per-group sum| = {s.abs().mean():.4f}   "
              f"(typical individual |asymmetry| = {sc:.4f})   "
              f"ratio = {s.abs().mean()/sc:.3f}")

    print("\n2. Is the asymmetry a stable individual trait? Intraclass correlation across days:")
    print("   (0 = pure day-to-day noise, 1 = perfectly stable individual property)\n")
    for f, r, name in PAIRS:
        a = icc(d, f"asym_{name}")
        fw = icc(d, f); rv = icc(d, r)
        print(f"   {name:14s} asymmetry ICC = {a:.3f}    "
              f"(forward alone {fw:.3f}, reverse alone {rv:.3f})")

    print("\n3. Hierarchy steepness per group, and its consistency across days:")
    for f, r, name in PAIRS[:2]:
        st = d.groupby(["GroupNumber", "Day"])[f"asym_{name}"].std()
        st = st.reset_index().rename(columns={f"asym_{name}": "steep"})
        i = icc(st, "steep", unit="GroupNumber")
        print(f"   {name:14s} mean steepness = {st.steep.mean():.4f}   "
              f"group-level ICC across days = {i:.3f}")
