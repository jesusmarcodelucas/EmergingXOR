#!/usr/bin/env python3
"""
ecl_sweep.py -- the two chloride set-point robustness sweeps for EmergingXOR.

Both set points of the two-phase model are favourable choices, and the appendix
now says so.  This closes the remaining [TO BE COMPLETED] by measuring how far
each can be moved before the result breaks.

  BUILDER sweep   ECl_immature swept, ECl_mature fixed at -65 mV.
                  The spike threshold is VT = -50 mV, so -40 and -45 mV are
                  frankly EXCITATORY while -55 and -60 mV are depolarising but
                  SUBTHRESHOLD.  Question: does the builder still lift the arms
                  over the firing cliff (w = 0.21) when GABA never drives the
                  laterals to threshold on its own?

  GOVERNOR sweep  ECl_mature swept, ECl_immature fixed at -40 mV.
                  -65 mV is exactly rest (purely shunting); -70 and -75 mV add a
                  subtractive, hyperpolarising component.  Question: does
                  coincidence suppression survive when inhibition is not purely
                  divisive?

phase_sim3.ECL_IMMATURE / ECL_MATURE are module-level constants that run_phases
resolves at call time, so they are patched per worker rather than by editing
phase_sim3.py.  Each run is classified INDIVIDUALLY and the distribution is
reported, never the classification of an average.

USAGE
    python3 ecl_sweep.py --seeds 5 --workers 32              # both sweeps
    python3 ecl_sweep.py --which builder --seeds 20 --workers 32
    python3 ecl_sweep.py --plot                              # figure from cache

Resumable: results are cached in ecl_sweep.json and completed jobs are skipped.
Runtime is roughly 2-3 min per run with the numpy target, so 5 values x 5 seeds
is about an hour on one core and a few minutes on 32.
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import argparse, json, os
from collections import Counter
from multiprocessing import Pool

import numpy as np

CACHE = "ecl_sweep.json"

# VT = -50 mV.  Below it the builder is depolarising but subthreshold.
BUILDER_ECL = [-60.0, -55.0, -50.0, -45.0, -40.0]     # mV, mature held at -65
GOVERNOR_ECL = [-75.0, -70.0, -65.0, -60.0]           # mV, immature held at -40

ECL_MATURE_REF = -65.0
ECL_IMMATURE_REF = -40.0
CLIFF = 0.21


def _one(job):
    """Run one (sweep, ECl, seed) point.  Imported inside the worker."""
    from brian2 import mV
    import two_phase_model as P

    which, ecl_imm, ecl_mat, seed_, tb, tr, tg, bal, sf = job
    P.ECL_IMMATURE = ecl_imm * mV          # patched before run_phases resolves them
    P.ECL_MATURE = ecl_mat * mV

    r = P.run_phases(mode="full", seed_=seed_, t_build=tb, t_ramp=tr, t_gov=tg,
                     balanced=bool(bal), hub_shared_frac=sf)
    key = f"{which}|{ecl_imm}|{ecl_mat}|{seed_}"
    return key, dict(which=which, ecl_imm=ecl_imm, ecl_mat=ecl_mat, seed=seed_,
                     cls=P.classify(r), arm=P.arm_verdict(r),
                     R10=float(r["R10"]), R01=float(r["R01"]),
                     R11=float(r["R11"]), R00=float(r["R00"]),
                     A10=float(r["A10"]), A01=float(r["A01"]),
                     xor_index=float(r["xor_index"]),
                     symmetry=float(r["symmetry"]),
                     w_mean=float(np.mean(r["w"])))


def save(res):
    with open(CACHE + ".tmp", "w") as fh:
        json.dump(res, fh, indent=1)
    os.replace(CACHE + ".tmp", CACHE)


def cells(res, which, values):
    """Aggregate the per-seed runs of one sweep into per-ECl summaries."""
    out = []
    for v in values:
        rs = [d for d in res.values()
              if d["which"] == which and
              (d["ecl_imm"] if which == "builder" else d["ecl_mat"]) == v]
        if not rs:
            continue
        cls = [d["cls"] for d in rs]
        arm = [d["arm"] for d in rs]
        xi = np.array([d["xor_index"] for d in rs])
        out.append(dict(
            ecl=v, n=len(rs),
            xor=np.mean([c == "XOR" for c in cls]),
            orf=np.mean([c == "OR" for c in cls]),
            col=np.mean([c == "collapse" for c in cls]),
            grown=np.mean([a == "both grown" for a in arm]),
            w=float(np.median([d["w_mean"] for d in rs])),
            xi=(float(np.percentile(xi, 50)), float(np.percentile(xi, 25)),
                float(np.percentile(xi, 75))),
            R10=float(np.median([d["R10"] for d in rs])),
            R11=float(np.median([d["R11"] for d in rs])),
            dom=Counter(cls).most_common(1)[0][0]))
    return out


def table(res, which, values, label, held):
    rows = cells(res, which, values)
    if not rows:
        return rows
    print(f"\n=== {label}  ({held}) ===")
    hdr = (f"{'ECl':>7} | {'n':>3} | {'XOR':>5} {'OR':>5} {'coll':>5} | "
           f"{'grown':>5} | {'w_med':>6} | {'xor_index (med [IQR])':>24} | "
           f"{'R10':>6} {'R11':>6} | outcome")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        flag = "" if r["w"] > CLIFF else "  <- below cliff"
        print(f"{r['ecl']:7.0f} | {r['n']:3d} | {r['xor']:5.2f} {r['orf']:5.2f} "
              f"{r['col']:5.2f} | {r['grown']:5.2f} | {r['w']:6.2f} | "
              f"{r['xi'][0]:+6.2f} [{r['xi'][1]:+5.2f}, {r['xi'][2]:+5.2f}]     | "
              f"{r['R10']:6.1f} {r['R11']:6.1f} | {r['dom']}{flag}")
    return rows


def plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

    B = cells(res, "builder", BUILDER_ECL)
    G = cells(res, "governor", GOVERNOR_ECL)
    if not B and not G:
        print("nothing cached to plot"); return

    fig, ax = plt.subplots(1, 2, figsize=(13.0, 4.6))
    for a, rows, ttl, xlab, vline, vtxt in [
            (ax[0], B, "A  builder set point  ($E_{Cl}^{\\rm mature}=-65$ mV)",
             "$E_{Cl}$ during the builder phase (mV)", -50.0, "spike threshold"),
            (ax[1], G, "B  governor set point  ($E_{Cl}^{\\rm builder}=-40$ mV)",
             "$E_{Cl}$ during the governor phase (mV)", -65.0, "rest (pure shunt)")]:
        if not rows:
            a.set_visible(False); continue
        x = [r["ecl"] for r in rows]
        a.plot(x, [r["xor"] for r in rows], "-o", color="#2e7d32", lw=2.6, ms=7,
               label="XOR fraction")
        a.plot(x, [r["grown"] for r in rows], "-s", color="#1f4e79", lw=2.2, ms=6,
               label="both-grown fraction")
        a.plot(x, [r["col"] for r in rows], "-^", color="#9e9e9e", lw=1.8, ms=6,
               label="collapse fraction")
        a.axvline(vline, ls=":", color="#c0392b", lw=1.4)
        a.text(vline, 1.04, vtxt, fontsize=8.5, color="#c0392b", ha="center")
        a.set_ylim(-0.04, 1.12); a.set_xlabel(xlab, fontsize=11)
        a.set_ylabel("fraction of runs", fontsize=11.5)
        a.set_title(ttl, weight="bold", loc="left", fontsize=11.5)
        a2 = a.twinx()
        a2.plot(x, [r["w"] for r in rows], "--d", color="#8a2be2", lw=1.8, ms=5)
        a2.axhline(CLIFF, ls="--", color="#8a2be2", lw=1.0, alpha=0.5)
        a2.set_ylabel("median arm weight $w$", color="#8a2be2", fontsize=11)
        a2.tick_params(axis="y", labelcolor="#8a2be2")
        a2.set_ylim(0, max(0.7, max(r["w"] for r in rows) * 1.25))
        a.legend(fontsize=9.5, frameon=False, loc="center left")
        for s in ("top",):
            a.spines[s].set_visible(False); a2.spines[s].set_visible(False)
    n = (B or G)[0]["n"]
    fig.suptitle(f"Robustness of the two-phase result to the chloride set points "
                 f"(full sequence, n={n} seeds/point; dashed purple: arm weight, "
                 f"line at the firing cliff)", fontsize=12, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig("ecl_sweep_en.png", dpi=160, bbox_inches="tight", facecolor="white")
    print("wrote ecl_sweep_en.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="both",
                    choices=["builder", "governor", "both"])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--seed0", type=int, default=11)
    ap.add_argument("--t_build", type=float, default=18.0)
    ap.add_argument("--t_ramp", type=float, default=4.0)
    ap.add_argument("--t_gov", type=float, default=12.0)
    ap.add_argument("--balanced", type=int, default=1)
    ap.add_argument("--shared", type=float, default=1.0)
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--plot", action="store_true",
                    help="only redraw the figure from the cache")
    a = ap.parse_args()

    res = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    if a.plot:
        table(res, "builder", BUILDER_ECL, "BUILDER sweep", "mature fixed at -65 mV")
        table(res, "governor", GOVERNOR_ECL, "GOVERNOR sweep", "builder fixed at -40 mV")
        plot(res); return

    seeds = list(range(a.seed0, a.seed0 + a.seeds))
    jobs = []
    if a.which in ("builder", "both"):
        for e in BUILDER_ECL:
            for s in seeds:
                jobs.append(("builder", e, ECL_MATURE_REF, s,
                             a.t_build, a.t_ramp, a.t_gov, a.balanced, a.shared))
    if a.which in ("governor", "both"):
        for e in GOVERNOR_ECL:
            for s in seeds:
                jobs.append(("governor", ECL_IMMATURE_REF, e, s,
                             a.t_build, a.t_ramp, a.t_gov, a.balanced, a.shared))
    todo = [j for j in jobs
            if f"{j[0]}|{j[1]}|{j[2]}|{j[3]}" not in res]
    print(f"{len(jobs)} runs total, {len(todo)} to do, {a.workers} workers", flush=True)

    if todo:
        with Pool(a.workers) as pool:
            for i, (k, v) in enumerate(pool.imap_unordered(_one, todo), 1):
                res[k] = v
                save(res)
                print(f"[{i:3d}/{len(todo)}] {k} -> {v['cls']:8s} | "
                      f"arms {v['arm']:12s} w={v['w_mean']:.2f} | "
                      f"{v['R10']:5.0f}/{v['R01']:5.0f}/{v['R11']:5.0f} | "
                      f"xi={v['xor_index']:+.2f}", flush=True)

    table(res, "builder", BUILDER_ECL, "BUILDER sweep", "mature fixed at -65 mV")
    table(res, "governor", GOVERNOR_ECL, "GOVERNOR sweep", "builder fixed at -40 mV")
    print("\nReading the builder sweep: the manuscript's claim is that the")
    print("depolarising phase lifts the arms over the cliff (w > 0.21). Any row")
    print("with a median w at or near 0.10 means the builder never taught, and")
    print("locates the subthreshold boundary of the mechanism.")
    plot(res)


if __name__ == "__main__":
    main()
