#!/usr/bin/env python3
"""
make_phase_figures.py  --  regenerate the three two-phase figures used in the
appendix (app:growth) of EmergingXOR.tex, with the _en filenames the .tex expects:

    phase_trace_en.png      (Fig. builder->switch->governor dynamics)
    phase_factorial_en.png  (Fig. 5-way factorial: each phase + plasticity necessary)
    window2d_en.png         (Fig. developmental window, t_build x mature shunt)

DATA IT READS
-------------
  phase_traces2.npz  : clean recorded traces of the FULL sequence.
                       Produced by run_phase_experiments.py (which runs phase_sim2.py
                       with record=True).  Looked up in ./ then /tmp/.
  phase_factorial.json : the 5 factorial conditions + arm-acquisition verdicts.
                       Also produced by run_phase_experiments.py.  If absent, the
                       hardened reference values below are used (the ones in the PDF).
  window2d.json      : the 2D sweep produced by window2d.py.

So the usual order is:
    python3 window2d.py --seeds 3 --workers 8     # -> window2d.json   (slow, ~1-2 h)
    python3 run_phase_experiments.py              # -> phase_factorial.json + phase_traces2.npz
    python3 make_phase_figures.py                 # -> the three *_en.png
"""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import Counter

plt.rcParams.update({"font.family": "DejaVu Sans"})
COL = {"XOR": "#2e7d32", "OR": "#e08a00", "collapse": "#9e9e9e"}


def _find(name):
    for d in (".", "data", "/tmp"):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


# ----------------------------------------------------------------------------
# (1) phase_trace_en.png   -- builder lifts the arms over the cliff, switch, XOR
# ----------------------------------------------------------------------------
def trace_fig(tb=18.0, tr=4.0, tg=12.0, cliff=0.21):
    p = _find("phase_traces2.npz")
    if p is None:
        print("[trace]  phase_traces2.npz not found -- run run_phase_experiments.py first. Skipping.")
        return
    d = np.load(p)
    tw, W, tl, ECL = d["t_w"], d["W"], d["t_l"], d["ECL"]

    def band(ax):
        ax.axvspan(0, tb, color="#fde9d0", alpha=0.6, zorder=0)
        ax.axvspan(tb, tb + tr, color="#e0e0e0", alpha=0.7, zorder=0)
        ax.axvspan(tb + tr, tb + tr + tg, color="#d8ecd8", alpha=0.6, zorder=0)

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(10.4, 6.6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1, 2.4], hspace=0.08))
    band(a0); a0.plot(tl, ECL, color="#a11a1a", lw=2.6)
    a0.axhline(-50, ls=":", color="#999", lw=1.2)
    a0.text(tb + tr + tg, -49, "threshold", fontsize=9.5, color="#999", va="bottom", ha="right")
    a0.set_ylabel("$E_{Cl}$ (mV)", fontsize=12.5); a0.set_ylim(-68, -36)
    for s in ("top", "right"): a0.spines[s].set_visible(False)
    a0.tick_params(labelsize=10.5)

    band(a1)
    a1.plot(tw, W[0], color="#1f4e79", lw=2.8, label="arm $I\\to L_i$")
    a1.plot(tw, W[1], color="#3a8ec0", lw=2.8, label="arm $P\\to L_p$")
    a1.axhline(cliff, ls="--", color="#c0392b", lw=2.0)
    # cliff label moved to the left margin: it used to sit in the same corner as
    # the legend and the GOVERNOR annotation, and the three overlapped.
    a1.text(0.4, cliff + 0.015, "firing cliff", fontsize=10.5,
            color="#c0392b", ha="left", va="bottom")
    a1.set_ylabel("arm weight  $w$", fontsize=12.5)
    a1.set_xlabel("time (s)", fontsize=12.5); a1.set_xlim(0, tb + tr + tg)
    a1.set_ylim(-0.06, max(0.72, float(np.max(W)) * 1.28))   # headroom for the
    #                                        legend, floor for the phase captions
    a1.legend(fontsize=11.5, frameon=False, loc="upper left",
              bbox_to_anchor=(0.01, 0.99))
    for s in ("top", "right"): a1.spines[s].set_visible(False)
    a1.tick_params(labelsize=10.5)

    ph = [("BUILDER", "depolarising GABA (+)\n$\\to$ GDPs lift the arms", tb / 2, "#8a4b00"),
          ("SWITCH", "$E_{Cl}\\downarrow$", tb + tr / 2, "#555"),
          ("GOVERNOR", "shunting GABA ($-$)\n$\\to$ mature XOR", tb + tr + tg / 2, "#2e6b2e")]
    for head, sub, xc, cc in ph:
        a0.text(xc, -39, head, ha="center", va="top", fontsize=10.5, weight="bold", color=cc)
        a1.text(xc, -0.055, sub, ha="center", va="bottom", fontsize=9, color=cc)
    fig.suptitle("How the switch matures XOR: the builder lifts the arms over the cliff; the governor then uses them",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("phase_trace_en.png", dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[trace]  wrote phase_trace_en.png")


# ----------------------------------------------------------------------------
# (2) phase_factorial_en.png  -- 5 conditions + arm-acquisition verdicts
# ----------------------------------------------------------------------------
# Reference values from the run shown in the PDF (t_build=18, seed 11).
# KEPT FOR COMPARISON ONLY -- never used as a fallback, see factorial_fig().
# Order of vals = read-out rate for patterns [10, 01, 11, 00].
_FAC_DEFAULT = [
    ("governor only",        [0, 0, 0, 0],       "collapse", "[0.10, 0.10]", "neither"),
    ("no hub",               [0, 0, 0, 0],       "collapse", "[0.10, 0.10]", "neither"),
    ("builder only",         [119, 132, 340, 0], "OR",       "[0.77, 0.75]", "both"),
    ("no plasticity ($\\eta{=}0$)", [0, 0, 0, 0], "collapse", "[0.10, 0.10]", "neither"),
    ("FULL sequence",        [36, 38, 3, 0],     "XOR",      "[0.46, 0.45]", "both"),
]
_ACQ = {"both": "both grown", "neither": "neither grown", "asymmetric": "asymmetric"}


def factorial_fig():
    p = _find("phase_factorial.json")
    if p is None:
        print("[factorial]  phase_factorial.json NOT FOUND.")
        print("[factorial]  Refusing to fall back on hard-coded numbers: a figure of")
        print("[factorial]  constants is indistinguishable from a figure of results.")
        print("[factorial]  Run:  python3 run_phase_experiments.py")
        return
    fac = [(r["label"], r["vals"], r["cls"], r["arms"], r["acq"])
           for r in json.load(open(p))]
    print("[factorial]  using phase_factorial.json")

    patlab = ["10", "01", "11", "00"]; barcol = ["#1f4e79", "#3a6ea5", "#a11a1a", "#bbb"]
    fig, axes = plt.subplots(1, len(fac), figsize=(14.5, 4.6), sharey=True)
    for ax, (lab, vals, verdict, arms, acq) in zip(axes, fac):
        ax.bar(patlab, vals, color=barcol, edgecolor="k", lw=0.7)
        for i, v in enumerate(vals):
            ax.text(i, v + 7, f"{v:.0f}", ha="center", fontsize=10.5, weight="bold")
        ax.set_title(lab, fontsize=12, weight="bold")
        # Display wording follows the classification of App. A.3 and
        # analysis/rescore_runs.py: the class is "failure to teach".  The JSON
        # key stays "collapse" because that is what the cached runs recorded.
        vt = {"collapse": "failure to teach", "OR": "OR", "XOR": "XOR"}[verdict]
        ax.text(0.5, 0.93, vt, transform=ax.transAxes, ha="center", va="top", fontsize=12,
                weight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.28", fc=COL[verdict], ec="none"))
        ax.text(0.5, -0.20, f"arms {arms}\narms: {_ACQ.get(acq, acq)}",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=9, color="#555")
        ax.set_ylim(0, 380)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        # The axis label has to clear the two-line "arms ..." annotation above,
        # which is drawn in axes coordinates and so does not reserve any space.
        # labelpad is in points and is measured from the tick labels.
        ax.set_xlabel("input I,P", fontsize=10.5, labelpad=42)
    axes[0].set_ylabel("read-out rate (Hz)", fontsize=11.5)
    fig.suptitle("Each phase \u2014 and plasticity \u2014 is necessary; silencing the hub confirms which arms grew\n"
                 "heterogeneous hub, clean reset before probing, w0 = 0.1\n"
                 "builder only: OR at immature E_Cl; the same trained arms give XOR when re-probed at mature E_Cl",
                 fontsize=12, weight="bold")
    fig.tight_layout(rect=[0, 0.10, 1, 0.9])
    fig.savefig("phase_factorial_en.png", dpi=155, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[factorial]  wrote phase_factorial_en.png")


# ----------------------------------------------------------------------------
# (3) window2d_en.png  -- developmental window (t_build x mature shunt)
# ----------------------------------------------------------------------------
def window_fig():
    """Retired.  The developmental-window figure is now produced by window2d.py,
    which re-scores the labels (gating on the firing cliff rather than on an
    absolute read-out rate) and knows about the extended builder-duration grid.
    Drawing it here would silently reproduce the superseded classification."""
    print("[window]  window2d_en.png is now produced by window2d.py:")
    print("[window]      python3 window2d.py --plot")
    print("[window]  skipping.")


if __name__ == "__main__":
    trace_fig()
    factorial_fig()
    window_fig()
