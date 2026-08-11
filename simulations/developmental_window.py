#!/usr/bin/env python3
"""
window2d.py -- the developmental window, re-done.

Changes from the version that produced the published Fig. 6:

  * imports phase_sim3 (the old file imported `phase_sim`, which no longer
    exists, so the script did not run at all);
  * seeds are generated rather than drawn from a hard-coded list of five, so
    --seeds 20 actually gives 20;
  * the cache stores w_mean and xor_index, so runs can be re-scored offline;
  * outcomes are labelled with rescore.py's threshold-free criteria (gating on
    whether the arms cleared the firing cliff) instead of an absolute read-out
    rate.  Re-scoring the old cache moved 22 of 26 "collapse" runs to XOR;
  * the abscissa is named honestly.  The hub->lateral weight enters the same
    conductance in both phases and only ECl changes between them, so it sets
    the strength of the depolarising builder as well as that of the mature
    shunt.  That is the thesis, not a bug, and the axis is labelled as the
    coupling rather than as "mature shunt strength".
  * --w_build runs the CONTROL in which the builder-phase weight is pinned
    while the governor-phase weight is swept, which separates the two roles.
    Leave it unset for the thesis condition.

    python3 window2d.py --seeds 20 --workers 8
    python3 window2d.py --seeds 20 --workers 8 --w_build 22   # control
    python3 window2d.py --plot
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import argparse, json, os
from collections import defaultdict

import numpy as np

T_BUILD = [8.0, 15.0, 20.0, 25.0, 32.0, 40.0, 50.0, 65.0]
# 50 and 65 s were added to test an extrapolation: within the original
# grid xor_index declines monotonically with builder duration but never
# reaches the OR criterion (minimum 0.60 against a threshold of 0.50).
# Column-wise linear fits put the crossing near 50-59 s.
W_HL = [10.0, 16.0, 22.0, 28.0, 34.0]
CLIFF = 0.21
XI_XOR = 0.5
# Probes last 1 s, so a rate in Hz is a spike count.  xor_index =
# (smin - R11)/smin is a ratio of small counts when smin is low, and its
# sampling error swamps the quantity: smin = 1 Hz with R11 = 0 gives 1.00 off
# two spikes.  MIN_COUNT does NOT classify anything -- classification is on the
# cliff, as before -- it only marks cells where the index is too imprecise to
# quote, so the figure shows a dash instead of a spuriously perfect value.
MIN_COUNT = 10.0


def cache_name(w_build):
    return "window2d.json" if w_build is None else f"window2d_wb{w_build:g}.json"


def one(job):
    from brian2 import nS
    import two_phase_model as P
    tb, whl, sd, w_build = job
    kw = {} if w_build is None else dict(w_hl_build=w_build * nS)
    r = P.run_phases(mode="full", t_build=tb, w_hl=whl * nS, seed_=sd, **kw)
    return (f"{tb}|{whl}|{sd}",
            dict(tb=tb, whl=whl, seed=sd,
                 cls_old=P.classify(r), arm=P.arm_verdict(r),
                 R10=r["R10"], R01=r["R01"], R11=r["R11"], R00=r["R00"],
                 A10=r["A10"], A01=r["A01"],
                 xor_index=float(r["xor_index"]),
                 w_mean=float(np.mean(r["w"])), w=[float(x) for x in r["w"]]))


def label(d):
    """Outcome, gated on the cliff and never on an absolute rate.

    The gate is applied per arm, not to their mean: with two arms a mean hides
    asymmetry (0.30 and 0.14 average to 0.22 and would clear a 0.21 cliff).
    """
    grown = [w > CLIFF for w in d["w"]]
    if not any(grown):
        return "not_taught"
    if not all(grown):
        return "asymmetric"
    if min(d["R10"], d["R01"]) <= 0:
        return "silent"                          # xor_index undefined here
    return "XOR" if d["xor_index"] >= XI_XOR else "OR"


ORDER = ("not_taught", "asymmetric", "silent", "OR", "XOR")   # fixed precedence for ties


def dominant(labs):
    """Majority class with a DETERMINISTIC tie-break, plus an explicit tie flag.

    max(set(labs), key=labs.count) breaks ties by set-iteration order, which
    depends on Python's string hash randomisation and therefore changes between
    runs on identical data.  Two cells of this sweep are exact 10/10 ties, and
    they flipped label between redraws before this was fixed.  Ties are now
    resolved by a fixed precedence AND reported, because an exact tie is a
    result about the cell, not a detail to hide.
    """
    counts = {l: labs.count(l) for l in set(labs)}
    top = max(counts.values())
    winners = [l for l in ORDER if counts.get(l, 0) == top]
    return winners[0], top, len(labs), len(winners) > 1


def xi_is_reliable(ds):
    """Is the median single-input rate high enough for xor_index to mean much?"""
    return float(np.median([min(d["R10"], d["R01"]) for d in ds])) >= MIN_COUNT


def save(res, path):
    with open(path + ".tmp", "w") as fh:
        json.dump(res, fh, indent=1)
    os.replace(path + ".tmp", path)


def report(res):
    cell = defaultdict(list)
    for d in res.values():
        cell[(d["tb"], d["whl"])].append(d)
    print("\n" + "=" * 78)
    print("outcome (re-scored) / median arm weight / median xor_index")
    print(f"{'t_build':>8} |" + "".join(f"{w:>13.0f} nS" for w in W_HL))
    print("-" * 78)
    for tb in sorted(T_BUILD, reverse=True):
        row = f"{tb:8.0f} |"
        for whl in W_HL:
            ds = cell.get((tb, whl), [])
            if not ds:
                row += f"{'-':>16}"; continue
            labs = [label(d) for d in ds]
            dom, cnt, tot, tied = dominant(labs)
            frac = cnt / tot
            w = np.median([d["w_mean"] for d in ds])
            xi = np.median([d["xor_index"] for d in ds])
            xs = f"{xi:+4.2f}" if xi_is_reliable(ds) else "  --"
            mark = "=" if tied else " "
            row += f"{dom[:3].upper():>4}{frac:4.2f}{mark}{w:4.2f}/{xs}"
        print(row)
    n_old = sum(1 for d in res.values() if d.get("cls_old") != label(d))
    print(f"\n{len(res)} runs; {n_old} would be labelled differently by the "
          f"old absolute-rate classifier")
    ties = []
    for (tb, whl), ds in sorted(cell.items()):
        if dominant([label(d) for d in ds])[3]:
            ties.append((tb, whl))
    if ties:
        print(f"exact ties (marked '=' above, TIE in the figure): {ties}")
        print("  those cells are genuinely undecided and should be reported as such.")

    # Long builders push the arms towards the hard clip at WMAX.  A cell whose
    # arms sit at the ceiling is not measuring over-growth any more, it is
    # measuring WMAX, so flag it rather than let it enter the figure silently.
    try:
        from phase_sim3 import WMAX
    except Exception:
        WMAX = 2.6
    near = [d for d in res.values() if d["w_mean"] >= 0.95 * WMAX]
    if near:
        cells = sorted({(d["tb"], d["whl"]) for d in near})
        print(f"WARNING: {len(near)} runs have arms within 5% of the clip "
              f"WMAX={WMAX}; affected cells: {cells}")
        print("         those cells report the clip, not the model.")


def plot(res, out="window2d_en.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    COL = {"XOR": "#2e7d32", "OR": "#e08a00", "asymmetric": "#7b5ea7",
           "not_taught": "#9e9e9e", "silent": "#616161"}
    cell = defaultdict(list)
    for d in res.values():
        cell[(d["tb"], d["whl"])].append(d)

    fig, ax = plt.subplots(figsize=(9.6, 7.0))
    XI = np.full((len(T_BUILD), len(W_HL)), np.nan)
    for iy, tb in enumerate(T_BUILD):
        for ix, whl in enumerate(W_HL):
            ds = cell.get((tb, whl), [])
            if not ds:
                continue
            labs = [label(d) for d in ds]
            dom, n, _tot, tied = dominant(labs)
            xi = float(np.median([d["xor_index"] for d in ds]))
            w = float(np.median([d["w_mean"] for d in ds]))
            smin = float(np.median([min(d["R10"], d["R01"]) for d in ds]))
            ok = xi_is_reliable(ds)
            XI[iy, ix] = xi if ok else np.nan
            ax.add_patch(plt.Rectangle((ix - .5, iy - .5), 1, 1,
                                       fc=COL[dom], ec="white", lw=2.5))
            ax.text(ix, iy + .20, dom, ha="center", va="center",
                    fontsize=11, weight="bold", color="white")
            ax.text(ix, iy - .02, f"{n}/{len(ds)}" + ("  TIE" if tied else ""),
                    ha="center", va="center", fontsize=9, color="white")
            xs = f"xi={xi:+.2f}" if ok else "xi n/a"
            ax.text(ix, iy - .26, f"w={w:.2f}  {xs}  {smin:.0f} Hz",
                    ha="center", va="center", fontsize=7.6, color="white")
    ax.set_xticks(range(len(W_HL)))
    ax.set_xticklabels([f"{w:.0f}" for w in W_HL], fontsize=12)
    ax.set_yticks(range(len(T_BUILD)))
    ax.set_yticklabels([f"{t:.0f}" for t in T_BUILD], fontsize=12)
    ax.set_xlim(-.5, len(W_HL) - .5); ax.set_ylim(-.5, len(T_BUILD) - .5)
    ax.set_xlabel("hub$\\to$lateral coupling (nS per synapse)\n"
                  "acts as depolarising drive in the builder phase and as "
                  "shunt in the governor phase", fontsize=11.5)
    ax.set_ylabel("builder duration $t_{\\mathrm{build}}$ (s) "
                  "\u2014 when the switch happens", fontsize=12)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.legend(handles=[Patch(fc=COL[k], label=k)
                       for k in ("XOR", "OR", "asymmetric", "not_taught")],
              loc="upper center", bbox_to_anchor=(.5, -.13), ncol=3,
              frameon=False, fontsize=11)
    fig.suptitle("A developmental window for the switch\n"
                 "(cell: re-scored class, vote, median arm weight and "
                 "xor_index)", fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, .94])
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    from multiprocessing import Pool
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--seed0", type=int, default=11)
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--w_build", type=float, default=None,
                    help="CONTROL: pin the builder-phase hub weight (nS)")
    ap.add_argument("--plot", action="store_true",
                    help="only redraw from the cache")
    a = ap.parse_args()

    path = cache_name(a.w_build)
    res = json.load(open(path)) if os.path.exists(path) else {}
    if a.plot:
        report(res); plot(res); raise SystemExit

    seeds = list(range(a.seed0, a.seed0 + a.seeds))
    jobs = [(tb, whl, sd, a.w_build) for tb in T_BUILD for whl in W_HL
            for sd in seeds if f"{tb}|{whl}|{sd}" not in res]
    total = len(T_BUILD) * len(W_HL) * len(seeds)
    print(f"{total} runs, {len(jobs)} to do, {a.workers} workers "
          f"| cache {path}"
          + ("" if a.w_build is None else
             f" | CONTROL w_build={a.w_build} nS"), flush=True)
    if jobs:
        with Pool(a.workers) as pool:
            for i, (k, v) in enumerate(pool.imap_unordered(one, jobs), 1):
                res[k] = v; save(res, path)
                print(f"[{i:4d}/{len(jobs)}] {k:>18s} -> {label(v):10s} "
                      f"w={v['w_mean']:.2f} xi={v['xor_index']:+.2f}",
                      flush=True)
    report(res); plot(res)
