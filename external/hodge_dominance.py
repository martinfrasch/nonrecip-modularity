"""Exact Hodge decomposition of a dominance flow: is the directed structure a gradient?

This is the social analogue of the central question in the colloid work, and on a finite graph it
can be answered EXACTLY rather than approximately.

Given a group of n animals and the net directed interaction counts f_ij = -f_ji (say, chases from
i to j minus chases from j to i), the flow is an antisymmetric edge function. On the complete
graph with all triangles filled, the discrete Hodge decomposition is

    f  =  grad(s)  +  curl-component ,        grad(s)_ij = s_j - s_i

with no harmonic part. The two pieces have direct interpretations:

  grad(s)  a TRANSITIVE hierarchy. The flow derives from a scalar rank potential s, exactly as a
           conservative force derives from a potential. If A beats B and B beats C then A beats C
           by construction.
  curl     an INTRANSITIVE, cyclic residual: rock-paper-scissors structure. This part has no
           potential and cannot be represented by any ranking.

The cyclic fraction ||curl|| / ||f|| is therefore a measured, unambiguous answer to "is this
antisymmetric coupling a gradient?" -- the question our manuscript notes is conflated with
antisymmetry itself (Sec. 1). In the colloid system we could only test the consequences of that
identification; here it can be computed directly.

The analysis needs the ORDERED pair counts, which the published aggregate table does not contain.
This module implements it and validates it on constructed cases, so the data request is for a
ready analysis rather than a proposal.
"""
from __future__ import annotations
import numpy as np


def hodge(f):
    """Decompose an antisymmetric flow matrix into gradient and cyclic parts.

    Returns (s, grad_part, curl_part, cyclic_fraction).
    s is the least-squares rank potential; grad_part_ij = s_j - s_i.
    """
    f = np.asarray(f, float)
    n = f.shape[0]
    f = 0.5 * (f - f.T)                     # enforce antisymmetry
    # Least-squares node potential minimising ||f - grad(s)||^2 with grad(s)_ij = s_j - s_i.
    # Setting the derivative to zero and fixing the gauge sum(s) = 0 gives
    #     s_k = -(1/n) sum_j f_kj
    # i.e. MINUS the row mean. (An earlier version dropped the sign, which returned a cyclic
    # fraction of 2.0 for a perfectly transitive flow and reversed the recovered ranks.)
    s = -f.mean(1)
    s = s - s.mean()
    g = s[None, :] - s[:, None]              # grad(s)_ij = s_j - s_i
    c = f - g
    nf = np.linalg.norm(f)
    return s, g, c, (np.linalg.norm(c) / nf if nf > 0 else 0.0)


def hodge_logodds(wins, prior=0.5):
    """Hodge decomposition of the LOG-ODDS flow, which is the correct scale.

    A Bradley-Terry hierarchy has logit(p_ij) = r_i - r_j exactly, so the log-odds flow is
    exactly a gradient. The raw net-count flow is not: the count difference is a logistic
    function of the rank difference, so a perfectly transitive Bradley-Terry group retains a
    spurious cyclic fraction of about 0.16 under a linear decomposition of counts, independent
    of how much data is collected. Decomposing log-odds removes that confound, so the residual
    cyclic fraction reflects genuine intransitivity rather than the nonlinearity of the
    win-probability function.
    """
    w = np.asarray(wins, float)
    n = w.shape[0]
    tot = w + w.T
    p = (w + prior) / (tot + 2 * prior)
    f = np.log(p / (1 - p))
    np.fill_diagonal(f, 0.0)
    f = np.where(tot > 0, f, 0.0)
    return hodge(f)


def _demo():
    print("=== validation on constructed cases (n = 4) ===\n")
    # perfectly transitive: rank 3 > 2 > 1 > 0, flow proportional to rank difference
    r = np.array([0., 1., 2., 3.])
    T = r[None, :] - r[:, None]
    s, g, c, cf = hodge(T)
    print(f"  perfectly transitive      cyclic fraction = {cf:.4f}   (expect 0)")
    print(f"     recovered ranks {np.round(s - s.min(), 3)}  (true {r})")

    # pure 4-cycle: 0->1->2->3->0, no ranking can represent it
    C = np.zeros((4, 4))
    for i in range(4):
        C[i, (i + 1) % 4] = 1.0
        C[(i + 1) % 4, i] = -1.0
    s, g, c, cf = hodge(C)
    print(f"\n  pure 4-cycle              cyclic fraction = {cf:.4f}   (expect 1)")
    print(f"     recovered ranks {np.round(s, 3)}  (flat, as it must be)")

    # mixture
    for w in (0.25, 0.5, 0.75):
        M = (1 - w) * T / np.linalg.norm(T) + w * C / np.linalg.norm(C)
        _, _, _, cf = hodge(M)
        print(f"\n  mixture, cyclic weight {w:.2f}   cyclic fraction = {cf:.4f}")

    # what a noisy 4-animal group might look like
    rng = np.random.default_rng(0)
    print("\n=== sampling noise floor: transitive truth + Poisson counting noise ===")
    print("   interactions/pair   cyclic frac (raw counts)   cyclic frac (log-odds)")
    for lam in (5, 20, 100, 500):
        cfs, cls = [], []
        for _ in range(200):
            p = 1 / (1 + np.exp(-(r[:, None] - r[None, :])))   # Bradley-Terry: i beats j
            w = rng.binomial(lam, np.clip(p, 0, 1))
            cfs.append(hodge(w - w.T)[3])
            cls.append(hodge_logodds(w)[3])
        print(f"        {lam:5d}         {np.mean(cfs):.3f} +- {np.std(cfs):.3f}"
              f"        {np.mean(cls):.3f} +- {np.std(cls):.3f}")
    print("\n   Raw counts retain a floor near 0.16 however much data is collected, because a")
    print("   Bradley-Terry hierarchy is logistic in the rank difference while the Hodge")
    print("   gradient is linear. Log-odds removes that confound and the floor falls with N,")
    print("   so log-odds is the correct scale and the per-pair interaction count must be")
    print("   reported to establish the residual noise floor.")


if __name__ == "__main__":
    _demo()
