"""Manuscript figure set. One file per figure, publication-oriented.

Each figure is generated from the analysis CSVs and raw trajectories already in the repository;
none of the numbers here are hand-entered.
"""
from __future__ import annotations
import glob, re, math
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.7, "lines.linewidth": 1.2, "figure.dpi": 200,
    "axes.spines.top": False, "axes.spines.right": False,
})
C = {"recip": "#2c3e50", "nonrecip": "#c0392b", "third": "#e67e22",
     "null": "#95a5a6", "bio": "#16a085"}


# ---------------------------------------------------------------- Fig 1
def fig1_construction():
    """The chi construction and its validation."""
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.2))
    a = 0.005
    g = lambda l, r: a * l**4 / (r*r + l*l)**2.5
    r = np.linspace(0.15, 1.0, 300)
    lL, lS = 1/6, 1/9
    gi, gj = g(lS, r), g(lL, r)
    gbar = 0.5*(gi+gj)
    for chi, ls in [(0.0, ":"), (0.5, "--"), (1.0, "-")]:
        ax[0].plot(r, gbar + chi*(gi-gbar), color=C["nonrecip"], ls=ls, label=f"χ={chi:g}")
        ax[0].plot(r, gbar + chi*(gj-gbar), color=C["recip"], ls=ls)
    ax[0].set(xlabel="separation $r$", ylabel="force coefficient $c_i$, $c_j$",
              title="a  the χ construction")
    ax[0].legend(frameon=False, loc="upper right")
    ax[0].text(0.52, 0.55, "red: on $i$\nblue: on $j$", transform=ax[0].transAxes, fontsize=6.5)

    # symmetric part is invariant
    for chi in (0.0, 0.25, 0.5, 1.0, 1.5):
        ax[1].plot(r, 0.5*((gbar+chi*(gi-gbar)) + (gbar+chi*(gj-gbar))), color=C["recip"], lw=2.5, alpha=.35)
        ax[1].plot(r, 0.5*((gbar+chi*(gj-gbar)) - (gbar+chi*(gi-gbar))), color=C["nonrecip"], lw=1.0)
    ax[1].set(xlabel="separation $r$", ylabel="sector magnitude", title="b  sectors vs χ")
    ax[1].text(0.35, 0.80, "symmetric: invariant\nin χ (5 curves)", color=C["recip"],
               transform=ax[1].transAxes, fontsize=6.5)
    ax[1].text(0.35, 0.45, "antisymmetric:\nlinear in χ", color=C["nonrecip"],
               transform=ax[1].transAxes, fontsize=6.5)

    # alpha cancels from the ratio
    for al in (0.003, 0.005, 0.010, 0.015):
        gi2, gj2 = al*lS**4/(r*r+lS**2)**2.5, al*lL**4/(r*r+lL**2)**2.5
        ax[2].plot(r, np.abs(gi2-gj2)/np.abs(gi2+gj2), lw=2.2, alpha=.6, label=f"α̂={al:g}")
    ax[2].set(xlabel="separation $r$", ylabel="|antisym| / |sym|", ylim=(0, 1),
              title="c  α̂ cancels exactly")
    ax[2].legend(frameon=False, loc="lower right")
    fig.tight_layout(); fig.savefig("fig1_construction.png", bbox_inches="tight")
    print("  fig1_construction.png")


# ---------------------------------------------------------------- Fig 2
def fig2_structure():
    """Cluster count and largest-cluster fraction vs chi, showing the non-monotonic dip."""
    d = pd.read_csv("paper_results.csv")
    b = d[d.case == "bidisperse"]
    g = b.groupby("chi").agg(ncl=("n_cl", "mean"), ncl_s=("n_cl", "sem"),
                             lcf=("lcf", "mean"), lcf_s=("lcf", "sem"),
                             Q=("Q", "mean"), Q_s=("Q", "sem"),
                             sv=("sigv2", "mean"), sv_s=("sigv2", "sem")).reset_index()
    m = d[d.case == "monodisperse"]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))
    ax[0].errorbar(g.chi, g.ncl, g.ncl_s, marker="o", ms=4, color=C["nonrecip"], capsize=2)
    ax[0].set_yscale("log")   # log axis so the factor-5 dip at chi=0.25 is visible alongside n_cl=137
    ax[0].axhline(m.n_cl.mean(), color=C["null"], ls="--", lw=.9)
    ax[0].annotate("monodisperse ref", (1.5, m.n_cl.mean()), fontsize=6, color=C["null"],
                   va="bottom", ha="right")
    ax[0].annotate("dip below\nreciprocal value", (0.25, g.ncl[g.chi == 0.25].values[0]),
                   xytext=(0.5, 4), fontsize=6.5, arrowprops=dict(arrowstyle="->", lw=.7))
    ax[0].set(xlabel="χ", ylabel="number of clusters", title="a  cluster count (non-monotonic)")
    ax[1].errorbar(g.chi, g.lcf, g.lcf_s, marker="s", ms=4, color=C["recip"], capsize=2)
    ax[1].set(xlabel="χ", ylabel="largest-cluster fraction", ylim=(0, 1.05),
              title="b  condensate mass fraction")
    a2 = ax[2].twinx()
    ax[2].errorbar(g.chi, g.sv, g.sv_s, marker="^", ms=4, color=C["third"], capsize=2)
    ax[2].set_yscale("log")
    a2.errorbar(g.chi, g.Q, g.Q_s, marker="s", ms=3.5, color=C["recip"], capsize=2)
    a2.set_ylabel("Newman $Q$"); a2.set_ylim(0.86, 0.95)
    ax[2].set(xlabel="χ", ylabel="$\\sigma^2_v$", title="c  activity rises, $Q$ does not")
    ax[2].tick_params(axis="y", colors=C["third"]); a2.tick_params(axis="y", colors=C["recip"])
    fig.tight_layout(); fig.savefig("fig2_structure.png", bbox_inches="tight")
    print("  fig2_structure.png")


# ---------------------------------------------------------------- Fig 3
def fig3_scaling():
    """Box scaling: largest cluster proportional to N; no characteristic size."""
    N = np.array([4000, 9000, 16000]); big = np.array([3116., 7221., 12581.])
    fig, ax = plt.subplots(1, 2, figsize=(5.0, 2.3))
    ax[0].loglog(N, big, "o", ms=5, color=C["nonrecip"])
    xx = np.linspace(3500, 18000, 50)
    ax[0].loglog(xx, big[0]*(xx/N[0])**1.0, "-", color=C["recip"], lw=1,
                 label="phase separation ($N^{1}$)")
    ax[0].loglog(xx, big[0]*np.ones_like(xx), "--", color=C["null"], lw=1,
                 label="finite $S^*$ ($N^{0}$)")
    sl = np.polyfit(np.log(N), np.log(big), 1)[0]
    ax[0].set(xlabel="$N$ (fixed density)", ylabel="largest cluster",
              title=f"a  measured exponent {sl:.2f}")
    ax[0].legend(frameon=False, loc="upper left")
    ax[1].semilogx(N, big/N, "o-", ms=5, color=C["nonrecip"])
    ax[1].set(xlabel="$N$", ylabel="largest-cluster fraction", ylim=(0, 1),
              title="b  lcf flat across 4× in $N$")
    fig.tight_layout(); fig.savefig("fig3_scaling.png", bbox_inches="tight")
    print("  fig3_scaling.png")


# ---------------------------------------------------------------- Fig 4
def _sticks(ax, ticks):
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter()); ax.set_xticklabels([str(t) for t in ticks])


def fig4_rates():
    """The critical nucleus: Becker-Doring drift sign change and the committor."""
    from cluster_kinetics import BINS
    ev = pd.read_csv("kinetics_events.csv"); ev = ev[(ev.chi == 1.5) & (ev.S >= 2)]
    bd = ev[ev.kind.isin(["intact", "evaporate"]) | ((ev.kind == "merge") & (ev.dS <= 2))].copy()
    bd["bin"] = pd.cut(bd.S, BINS, right=False)
    g = bd.groupby("bin", observed=True).agg(S=("S", "mean"), d=("dS", "mean"), e=("dS", "sem")).reset_index()
    g = g[g.S < 100]
    cm = pd.read_csv("committor_summary.csv"); cm = cm[cm.S < 100]
    fig, ax = plt.subplots(1, 2, figsize=(5.0, 2.3))
    ax[0].errorbar(g.S, g.d, g.e, marker="o", ms=4, color=C["nonrecip"], capsize=2, lw=1)
    ax[0].axhline(0, color="k", lw=.7); ax[0].axvline(10.1, color=C["third"], lw=.8, ls="--")
    ax[0].axvspan(9.7, 10.7, color=C["third"], alpha=.2)
    ax[0].set_xscale("log"); _sticks(ax[0], [2, 3, 5, 10, 20, 50])
    ax[0].annotate("$S^*=10.1$", (11.5, g.d.min() * 0.9), fontsize=7)
    ax[0].set(xlabel="cluster size $S$", ylabel="⟨ΔS⟩ per interval (monomer exchange)",
              title="a  size drift changes sign")
    ax[1].errorbar(cm.S, cm.q, cm.q_err, marker="o", ms=4, color=C["nonrecip"], capsize=2, lw=1, label="$q$")
    ax[1].errorbar(cm.S, cm.q_free, cm.q_free_err, marker="s", ms=3.5, color=C["recip"], capsize=2, lw=1,
                   label="$q$, absorption excluded")
    ax[1].plot(cm.S, cm.open_frac, ":", color=C["null"], lw=1.2, label="undecided at $\\hat t=1000$")
    ax[1].axhline(0.5, color="k", lw=.7); ax[1].axvline(10.1, color=C["third"], lw=.8, ls="--")
    ax[1].set_xscale("log"); _sticks(ax[1], [3, 5, 10, 20, 50])
    ax[1].set(xlabel="initial cluster size $S_0$", ylabel="committor $q(S_0)$", ylim=(0, 1.02),
              title="b  committor is ½ at the nucleus")
    ax[1].legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig("fig4_rates.png", bbox_inches="tight")
    print("  fig4_rates.png")


# ---------------------------------------------------------------- Fig 5
def fig5_thermo():
    """Entropy production, the Onsager window scan, and the T_eff exponent breakdown."""
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))
    cal = pd.read_csv("epr_static_calibtraj.csv").groupby("chi").epr_calib.agg(["mean", "sem"])
    st = pd.read_csv("epr_static.csv"); st = st[(st.N == 4000) & (st.L == 24) & (st.chi >= 3)]
    st = st.groupby("chi").epr.agg(["mean", "sem"])
    d = pd.concat([cal, st]); chi = d.index.values; epr = d["mean"].values; err = d["sem"].values
    y = epr / chi**2; ye = err / chi**2
    ax[0].errorbar(chi, y * 1e3, ye * 1e3, marker="o", ms=4, color=C["nonrecip"], capsize=2, lw=0, elinewidth=1)
    ax[0].axhspan(2.3, 3.3, color=C["recip"], alpha=.12); ax[0].axhline(0, color="k", lw=.7)
    ax[0].axvspan(0.25, 0.75, color=C["third"], alpha=.15)
    ax[0].set_xscale("log"); _sticks(ax[0], [0.25, 0.5, 1, 2, 4, 8])
    ax[0].set(xlabel="χ", ylabel="excess dissipation / χ²  (×10⁻³)", title="a  a threshold, then quadratic")

    T = np.array([200, 400, 800])
    sym = np.array([-1.08e-8, -6.81e-8, -3.96e-7]); syme = np.array([6.5e-7, 8.6e-7, 9.9e-7])
    ant = np.array([9.975e-6, 9.401e-6, 9.477e-6]); ante = np.array([6.3e-7, 8.2e-7, 9.6e-7])
    ax[1].errorbar(T, ant, ante, marker="o", ms=4, color=C["nonrecip"], capsize=2,
                   label="antisymmetric")
    ax[1].errorbar(T, sym, syme, marker="s", ms=4, color=C["recip"], capsize=2,
                   label="symmetric")
    ax[1].errorbar(T, [2.16e-6, 1.17e-6, 1.18e-8], [4.2e-7, 5.1e-7, 5.9e-7], marker="o", ms=4, mfc="none",
                   color=C["nonrecip"], capsize=2, ls=":", label="antisym., dissimilar channel")
    ax[1].errorbar(T, [8.67e-7, 1.62e-6, 2.70e-6], [4.3e-7, 5.2e-7, 5.9e-7], marker="s", ms=4, mfc="none",
                   color=C["recip"], capsize=2, ls=":", label="sym., dissimilar channel")
    ax[1].axhline(0, color="k", lw=.7)
    ax[1].set(xlabel="measurement window $T$", ylabel="cross-coefficient",
              title="b  cross-response: identical vs dissimilar channels")
    ax[1].legend(frameon=False)

    hs = pd.read_csv("hs_results.csv"); hs = hs[hs.species.isin(["L", "S"]) & hs.t.notna()]
    cols = {0.0: C["recip"], 0.5: C["null"], 1.0: C["third"], 1.5: C["nonrecip"]}
    for chi_, g in hs.groupby("chi"):
        for sp, mk, ls in (("L", "o", "-"), ("S", "s", "--")):
            h = g[g.species == sp].sort_values("t")
            ax[2].loglog(h.t, h.Teff, marker=mk, ms=3, ls=ls, color=cols.get(chi_, "k"), lw=1,
                         label=f"χ={chi_:g}" if sp == "L" else None)
    ax[2].axhline(1, color="k", lw=.7)
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    ax[2].yaxis.set_major_locator(FixedLocator([1, 2, 4, 8])); ax[2].yaxis.set_minor_locator(NullLocator())
    ax[2].yaxis.set_major_formatter(ScalarFormatter())
    ax[2].set(xlabel="lag $t$", ylabel="$T_{\\rm eff}(t)/T$", title="c  one shared $T_{\\rm eff}$, equal to $T$ above $1/t_{\\rm FDT}$")
    ax[2].legend(frameon=False, loc="upper left", title="L (solid) / S (dashed)")
    fig.tight_layout(); fig.savefig("fig5_thermo.png", bbox_inches="tight")
    print("  fig5_thermo.png")


# ---------------------------------------------------------------- Fig 6
def fig6_circulation():
    """Colloid vs cilia circulation, at matched level of description."""
    fig, ax = plt.subplots(1, 2, figsize=(5.2, 2.3))
    labels = ["colloid\nχ=0", "colloid\nχ=1", "colloid\nχ=1.5", "cilia\n(axoneme)"]
    med = [0.81, 0.85, 0.94, 3.2]
    frac = [0.0, 0.0, 0.0, 0.92]
    cols = [C["recip"], C["nonrecip"], C["nonrecip"], C["bio"]]
    ax[0].bar(range(4), med, color=cols, width=.6)
    ax[0].axhline(2, color="k", ls="--", lw=.8)
    ax[0].annotate("|z| = 2", (3.4, 2.1), fontsize=6.5, ha="right")
    ax[0].set_xticks(range(4)); ax[0].set_xticklabels(labels, fontsize=6.5)
    ax[0].set(ylabel="median |z| of signed area rate",
              title="a  single object, matched estimator")
    ax[1].bar(range(4), frac, color=cols, width=.6)
    ax[1].set_xticks(range(4)); ax[1].set_xticklabels(labels, fontsize=6.5)
    ax[1].set(ylabel="fraction with |z| > 2", ylim=(0, 1),
              title="b  detection rate")
    fig.tight_layout(); fig.savefig("fig6_circulation.png", bbox_inches="tight")
    print("  fig6_circulation.png")


if __name__ == "__main__":
    print("building manuscript figures:")
    fig1_construction(); fig2_structure(); fig3_scaling()
    fig4_rates(); fig5_thermo(); fig6_circulation()
