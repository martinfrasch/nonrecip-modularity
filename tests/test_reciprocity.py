"""Check the chi knob does what the experiment design claims.

chi=0 must make the EHD force obey Newton's third law exactly (net pair force = 0,
so an isolated noise-free pair has a fixed centre of mass). chi=1 must reproduce
the original Hara et al. force. Run: python tests/test_reciprocity.py
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sim.simulate import _run


def pair_com_drift(chi, nsteps=2000, dt=0.05, L=6.0):
    """Two particles, unequal EHD radii, equal steric radii, zero noise."""
    pos = np.array([[3.0, 3.0], [3.25, 3.0]])
    svec = np.array([1 / 6, 1 / 6])          # equal steric -> equal mobility
    lvec = np.array([1 / 6, 1 / 9])          # unequal EHD  -> nonreciprocal at chi>0
    snaps = np.zeros((2, 2, 2)); snaps[0] = pos
    heat = np.zeros(2)
    com0 = pos.mean(axis=0).copy()
    _run(pos, svec, lvec, lvec**2, lvec**4, nsteps, dt, L, 0.005, chi,
         0.0, 1, nsteps, snaps, heat)          # sigma=0
    d = pos.mean(axis=0) - com0
    d -= L * np.round(d / L)
    return float(np.hypot(*d))


def main():
    d0 = pair_com_drift(0.0)
    d1 = pair_com_drift(1.0)
    dh = pair_com_drift(0.5)
    print(f"COM drift  chi=0.0 : {d0:.3e}   (must be ~0: reciprocal)")
    print(f"COM drift  chi=0.5 : {dh:.3e}")
    print(f"COM drift  chi=1.0 : {d1:.3e}   (must be >0: nonreciprocal)")

    assert d0 < 1e-12, f"chi=0 is not reciprocal: COM drifted {d0:.3e}"
    assert d1 > 1e-4, f"chi=1 shows no self-propulsion: {d1:.3e}"
    # nonreciprocal drive is linear in chi -> so is the drift, at fixed geometry
    assert abs(dh / d1 - 0.5) < 0.15, f"drift not ~linear in chi: {dh/d1:.3f} vs 0.5"

    # symmetric part must not depend on chi: identical particles => identical dynamics
    def same_particles(chi):
        pos = np.array([[3.0, 3.0], [3.3, 3.0]])
        s = np.array([1 / 6, 1 / 6]); l = np.array([1 / 6, 1 / 6])
        snaps = np.zeros((2, 2, 2)); snaps[0] = pos
        heat = np.zeros(2)
        _run(pos, s, l, l**2, l**4, 500, 0.05, 6.0, 0.005, chi, 0.0, 1, 500, snaps, heat)
        return pos.copy()
    a, b = same_particles(0.0), same_particles(1.0)
    assert np.allclose(a, b, atol=1e-14), "chi changed dynamics of an identical pair"
    print("equal-radius pair: chi=0 and chi=1 trajectories identical  (symmetric part invariant)")
    print("\nOK")


if __name__ == "__main__":
    main()
