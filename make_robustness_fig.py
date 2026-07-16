#!/usr/bin/env python3
"""
make_robustness_fig.py -- build the hub-independence sweep figure from the
robust_shared*.json files written by run_robustness.py.

It reads every robust_shared*.json in the current directory, keeps the
BALANCED runs (balanced==1) for the main sweep, sorts them by hub-sharing
fraction, and plots three panels vs shared_frac:

  A  XOR fraction  and  both-grown fraction        (the discrete outcome)
  B  xor_index median with IQR band                (the continuous selectivity)
  C  coincidence rate R11 median with IQR band     (what leaks in as the hub
                                                     de-correlates)

Output: robustness_sweep_en.png

Any fully-random (balanced==0) point is drawn as a separate reference marker
at its shared_frac.  Run it after the sweep finishes; it works with whatever
JSONs are present.
"""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_XOR, C_GROW, C_IDX, C_R11 = "#2e7d32", "#1f4e79", "#8a2be2", "#c0392b"


def load():
    bal, rand = [], []
    for fn in glob.glob("data/robust*shared*.json") + glob.glob("robust*shared*.json"):
        try:
            S = json.load(open(fn))["summary"]
        except Exception as e:
            print("skip", fn, e); continue
        rec = dict(sf=float(S["shared_frac"]),
                   xor=float(S["xor_fraction"]),
                   grown=float(S["both_grown_fraction"]),
                   idx=S["xor_index"], r11=S["R11"], n=int(S["n"]),
                   bal=int(S.get("balanced", 1)))
        (bal if rec["bal"] == 1 else rand).append(rec)
    bal.sort(key=lambda r: r["sf"], reverse=True)
    rand.sort(key=lambda r: r["sf"], reverse=True)
    return bal, rand


def main():
    bal, rand = load()
    if not bal:
        print("No robust_shared*.json (balanced) found in this directory.")
        return
    sf = [r["sf"] for r in bal]
    xor = [r["xor"] for r in bal]
    grown = [r["grown"] for r in bal]
    idx_m = [r["idx"][0] for r in bal]; idx_lo = [r["idx"][1] for r in bal]; idx_hi = [r["idx"][2] for r in bal]
    r11_m = [r["r11"][0] for r in bal]; r11_lo = [r["r11"][1] for r in bal]; r11_hi = [r["r11"][2] for r in bal]
    n = bal[0]["n"]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.4))

    # A: discrete outcome fractions
    ax[0].plot(sf, xor, "-o", color=C_XOR, lw=2.6, ms=7, label="XOR fraction")
    ax[0].plot(sf, grown, "-s", color=C_GROW, lw=2.2, ms=6, label="both-grown fraction")
    ax[0].set_ylim(-0.03, 1.05); ax[0].set_ylabel("fraction of runs", fontsize=12)
    ax[0].set_title("A  outcome vs hub sharing", weight="bold", loc="left", fontsize=12.5)
    ax[0].legend(fontsize=10.5, frameon=False, loc="lower right")

    # B: continuous selectivity index with IQR band
    ax[1].fill_between(sf, idx_lo, idx_hi, color=C_IDX, alpha=0.18)
    ax[1].plot(sf, idx_m, "-o", color=C_IDX, lw=2.6, ms=7)
    ax[1].axhline(0, ls="--", color="#888", lw=1.3)
    ax[1].text(min(sf), 0.02, "XOR", fontsize=9, color="#2e7d32", va="bottom")
    ax[1].text(min(sf), -0.05, "OR", fontsize=9, color="#e08a00", va="top")
    ax[1].set_ylabel("xor_index  (median, IQR)", fontsize=12)
    ax[1].set_title("B  selectivity vs hub sharing", weight="bold", loc="left", fontsize=12.5)

    # C: coincidence leak
    ax[2].fill_between(sf, r11_lo, r11_hi, color=C_R11, alpha=0.16)
    ax[2].plot(sf, r11_m, "-o", color=C_R11, lw=2.6, ms=7)
    ax[2].set_ylabel("coincidence rate $R_{11}$ (Hz)", fontsize=12)
    ax[2].set_title("C  coincidence leak vs hub sharing", weight="bold", loc="left", fontsize=12.5)

    # fully-random reference points (balanced==0), if any
    for r in rand:
        ax[0].plot(r["sf"], r["xor"], "^", color=C_XOR, ms=9, mfc="white", mew=1.6,
                   label="XOR (random exposure)" if r is rand[0] else None)
        ax[1].plot(r["sf"], r["idx"][0], "^", color=C_IDX, ms=9, mfc="white", mew=1.6)
    if rand:
        ax[0].legend(fontsize=9.5, frameon=False, loc="lower right")

    for a in ax:
        a.set_xlabel("hub input sharing  (shared_frac)\n1.0 = fully shared   \u2192   0.5 = partly independent", fontsize=11)
        a.set_xlim(max(sf) + 0.03, min(sf) - 0.03)   # 1.0 on the left
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(labelsize=10.5)

    fig.suptitle(f"XOR needs a correlated hub: coincidence suppression degrades as hub afferents become independent  "
                 f"(full sequence, n={n} seeds/point)", fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("robustness_sweep_en.png", dpi=160, bbox_inches="tight", facecolor="white")
    print("wrote robustness_sweep_en.png")

    print("\nshared_frac | XOR | both_grown | xor_index(med) | R11(med)")
    for r in bal:
        print(f"   {r['sf']:.2f}    | {r['xor']:.2f}|   {r['grown']:.2f}    |     {r['idx'][0]:+.2f}     |  {r['r11'][0]:.0f}")


if __name__ == "__main__":
    main()
