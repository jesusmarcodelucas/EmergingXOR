#!/usr/bin/env python3
"""
make_ecl_fig.py -- the ECl set-point figure, drawn from RE-SCORED labels.

The figure ecl_sweep.py produces natively uses phase_sim3.classify(), whose
absolute 25 Hz rate gate mislabels every successful run at ECl_build = -45 mV as
"collapse" (see rescore.py).  Drawn that way the XOR curve sits at zero exactly
where the circuit is most selective, which is the opposite of the truth.  This
script re-scores first and then plots.

Panel A (builder)  : outcome vs ECl_build, with median arm weight and the firing
                     cliff, so the reader can see that below VT the arms never
                     move at all.
Panel B (governor) : outcome vs ECl_mature, with the coincidence rate R11, which
                     is the mechanistically interesting quantity: it falls
                     monotonically as inhibition becomes more hyperpolarising.

    python3 make_ecl_fig.py                    # -> ecl_sweep_en.png
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import argparse, json
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rescore_runs import rescore, CLIFF

BUILDER_ECL = [-60.0, -55.0, -50.0, -45.0, -40.0]
GOVERNOR_ECL = [-75.0, -70.0, -65.0, -60.0]
VT = -50.0
EL = -65.0

C_XOR, C_GROWN, C_NT = "#2e7d32", "#1f4e79", "#9e9e9e"
C_W, C_R11 = "#8a2be2", "#c0392b"


def load(path):
    raw = json.load(open(path))
    grp = defaultdict(list)
    for d in raw.values():
        new, _ = rescore(d)
        d = dict(d, new=new)
        v = d["ecl_imm"] if d["which"] == "builder" else d["ecl_mat"]
        grp[(d["which"], v)].append(d)
    return grp


def stats(rows):
    n = len(rows)
    f = lambda lab: sum(r["new"] == lab for r in rows) / n
    return dict(n=n, xor=f("XOR"), orr=f("OR"), nt=f("not_taught"),
                w=float(np.median([r["w_mean"] for r in rows])),
                xi=[float(np.percentile([r["xor_index"] for r in rows], q))
                    for q in (50, 25, 75)],
                R10=float(np.median([r["R10"] for r in rows])),
                R11=float(np.median([r["R11"] for r in rows])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="ecl_sweep.json")
    ap.add_argument("--out", default="ecl_sweep_en.png")
    a = ap.parse_args()

    grp = load(a.file)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, ax = plt.subplots(1, 2, figsize=(14.6, 5.0))

    # ---------------------------------------------------------------- builder
    xs = [v for v in BUILDER_ECL if ("builder", v) in grp]
    S = [stats(grp[("builder", v)]) for v in xs]
    a0 = ax[0]
    a0.plot(xs, [s["xor"] for s in S], "-o", color=C_XOR, lw=2.8, ms=8,
            label="XOR fraction (re-scored)", zorder=4)
    a0.plot(xs, [s["nt"] for s in S], "-^", color=C_NT, lw=2.0, ms=7,
            label="arms never taught", zorder=3)
    a0.axvline(VT, ls=":", color=C_R11, lw=1.8)
    a0.text(VT - 0.3, 1.02, "spike threshold $V_T$", fontsize=9.5,
            color=C_R11, ha="right")
    a0.set_ylim(-0.05, 1.10); a0.set_ylabel("fraction of runs", fontsize=12)
    a0.set_xlabel("$E_{Cl}$ during the builder phase (mV)", fontsize=12)
    a0.set_title("A  builder set point  ($E_{Cl}^{\\mathrm{mature}}=-65$ mV)",
                 loc="left", weight="bold", fontsize=12.5)
    a0.legend(fontsize=10, frameon=False, loc="center left")
    b0 = a0.twinx()
    b0.plot(xs, [s["w"] for s in S], "--D", color=C_W, lw=2.0, ms=6)
    b0.axhline(CLIFF, ls="--", color=C_W, lw=1.0, alpha=0.5)
    b0.text(xs[0], CLIFF + 0.012, "firing cliff", fontsize=9, color=C_W)
    b0.set_ylabel("median arm weight $w$", color=C_W, fontsize=12)
    b0.tick_params(axis="y", labelcolor=C_W); b0.set_ylim(0, 0.70)

    # --------------------------------------------------------------- governor
    xs2 = [v for v in GOVERNOR_ECL if ("governor", v) in grp]
    S2 = [stats(grp[("governor", v)]) for v in xs2]
    a1 = ax[1]
    a1.plot(xs2, [s["xor"] for s in S2], "-o", color=C_XOR, lw=2.8, ms=8,
            label="XOR fraction (re-scored)", zorder=4)
    xim = [s["xi"][0] for s in S2]
    a1.plot(xs2, xim, "-d", color=C_W, lw=2.2, ms=7, label="xor_index (median)")
    a1.fill_between(xs2, [s["xi"][1] for s in S2], [s["xi"][2] for s in S2],
                    color=C_W, alpha=0.15)
    a1.axvline(EL, ls=":", color=C_R11, lw=1.8)
    a1.text(EL - 0.15, 1.02, "rest $E_L$ (pure shunt)", fontsize=9.5,
            color=C_R11, ha="right")
    a1.set_ylim(-0.05, 1.10); a1.set_ylabel("fraction / index", fontsize=12)
    a1.set_xlabel("$E_{Cl}$ during the governor phase (mV)", fontsize=12)
    a1.set_title("B  governor set point  ($E_{Cl}^{\\mathrm{builder}}=-40$ mV)",
                 loc="left", weight="bold", fontsize=12.5)
    a1.legend(fontsize=10, frameon=False, loc="center left")
    b1 = a1.twinx()
    b1.plot(xs2, [s["R11"] for s in S2], "-s", color=C_R11, lw=2.2, ms=7)
    b1.set_ylabel("coincidence rate $R_{11}$ (Hz)", color=C_R11, fontsize=12)
    b1.tick_params(axis="y", labelcolor=C_R11); b1.set_ylim(0, 16)

    for axx in (a0, a1):
        for s in ("top",):
            axx.spines[s].set_visible(False)

    fig.suptitle("Chloride set-point robustness of the two-phase model "
                 "(full sequence, $n=20$ seeds per point; labels re-scored on "
                 "the firing cliff, not on an absolute rate)",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(a.out, dpi=160, bbox_inches="tight", facecolor="white")
    print("wrote", a.out)

    print("\nbuilder:")
    for v, s in zip(xs, S):
        print(f"  {v:6.0f} mV | XOR {s['xor']:.2f} | not_taught {s['nt']:.2f} "
              f"| w {s['w']:.2f} | xi {s['xi'][0]:+.2f} | R10 {s['R10']:.0f} "
              f"| R11 {s['R11']:.0f}")
    print("governor:")
    for v, s in zip(xs2, S2):
        print(f"  {v:6.0f} mV | XOR {s['xor']:.2f} | w {s['w']:.2f} "
              f"| xi {s['xi'][0]:+.2f} [{s['xi'][1]:+.2f}, {s['xi'][2]:+.2f}] "
              f"| R10 {s['R10']:.0f} | R11 {s['R11']:.1f}")


if __name__ == "__main__":
    main()
