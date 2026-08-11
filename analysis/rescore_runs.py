#!/usr/bin/env python3
"""
rescore.py -- re-label completed runs from cached JSON, without re-simulating.

WHY
---
phase_sim3.classify() gates on an ABSOLUTE read-out rate:

    if min(R10, R01) < theta_R:   # theta_R = 25 Hz
        return "collapse"

That conflates two different failures:

  (a) the builder never taught          -> arms stay at w0, read-out silent
  (b) the builder taught, but the grown circuit operates at a LOWER RATE
      than 25 Hz while still suppressing coincidence cleanly

Case (b) is mislabelled "collapse" even when xor_index = 1.00.  It showed up
first in the ECl sweep: at ECl_build = -45 mV all 20 seeds grew their arms to
w ~ 0.31 (cliff = 0.21) and reached xor_index median 0.86 -- better selectivity
than the -40 mV headline condition -- yet every run was labelled collapse
because min(R10,R01) sat at 15-24 Hz.

arm_verdict() is NOT a safe substitute: it thresholds A10/A01 at the same
absolute 25 Hz, and at least one seed reports "neither grown" with w = 0.31.

THE FIX
-------
Gate on quantities that are threshold-free or circuit-intrinsic:

  taught     w_mean > CLIFF (0.21)      did the builder move the arms at all?
  active     min(R10,R01) > 0           does the read-out respond to singles?
  selective  xor_index >= XI_XOR (0.5)  same selectivity criterion as before
                                        (the original alpha=0.5 is exactly
                                         xor_index >= 0.5), now applied only
                                         once the circuit is known to be taught

Classes: not_taught / silent / OR / XOR, plus a leaky flag when R00 > 10 Hz.
Rate is reported as a DESCRIPTOR, never as a gate.

USAGE
    python3 rescore.py                       # ecl_sweep.json
    python3 rescore.py --file window2d.json  # check the Fig. 6 sweep too
    python3 rescore.py --file ecl_sweep.json --csv rescored.csv
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import argparse, json, os
from collections import Counter, defaultdict

import numpy as np

CLIFF = 0.21
XI_XOR = 0.5
LEAK = 10.0
OLD_THETA_R = 25.0


def rescore(d):
    """Return (new_class, leaky_flag) for one cached run dict.

    The cliff gate is applied to the INDIVIDUAL arm weights, not to their mean.
    With only two arms a mean or median hides asymmetry: w_I = 0.30 with
    w_P = 0.14 averages to 0.22 and would clear a 0.21 cliff although one arm
    never grew.  Runs where exactly one arm cleared the cliff get their own
    class rather than being absorbed into either success or failure.
    """
    def num(key, default=None):
        """Some older caches store null for fields they could not serialise, so
        a plain .get(key, default) still returns None.  Coerce explicitly."""
        v = d.get(key)
        return default if v is None else float(v)

    R10, R01, R11 = num("R10", 0.0), num("R01", 0.0), num("R11", 0.0)
    R00 = num("R00", 0.0)
    smin = min(R10, R01)
    leaky = R00 > LEAK

    # Whether each arm was taught, by the best evidence the file carries:
    #   1. the two arm weights against the cliff (best);
    #   2. the arm-acquisition probes A10/A01, measured with the hub silenced,
    #      where a non-zero response means that arm alone drives its lateral --
    #      which is what the cliff is defined as, so this is a per-arm test and
    #      not an arbitrary rate gate;
    #   3. the mean weight, which cannot see asymmetry.
    arms = d.get("w")
    if arms is not None:
        grown = [w > CLIFF for w in arms]
    elif num("A10") is not None and num("A01") is not None:
        grown = [num("A10") > 0.0, num("A01") > 0.0]
    elif num("w_mean") is not None:
        grown = [num("w_mean") > CLIFF] * 2
    else:
        # Nothing in the file records whether the arms were taught.  Rather than
        # return "unknown" for every run, fall back to selectivity alone and say
        # so loudly: for a sweep whose companion text already reports both arms
        # grown in every run, this still answers the question the re-scoring was
        # meant to answer, namely how many runs the old 25 Hz gate discarded.
        grown = [True, True]
    if not any(grown):
        return "not_taught", leaky
    if not all(grown):
        return "asymmetric", leaky
    if smin <= 0:
        return "silent", leaky                   # xor_index undefined here
    xi = num("xor_index")
    if xi is None:
        xi = (smin - R11) / smin
    return ("XOR" if xi >= XI_XOR else "OR"), leaky


def group_key(k, d):
    """Best-effort grouping key for whichever sweep produced the file."""
    if "which" in d:                       # ecl_sweep.json
        v = d["ecl_imm"] if d["which"] == "builder" else d["ecl_mat"]
        return (d["which"], v)
    if "ecl" in d:
        return (d.get("axis", "ecl"), d["ecl"])
    # window2d.json keys are "t_build|w_hl|seed" with no descriptive fields
    parts = str(k).split("|")
    if len(parts) == 3:
        try:
            return (f"tb={float(parts[0]):g}", float(parts[1]))
        except ValueError:
            pass
    return ("all", None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="ecl_sweep.json")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    if not os.path.exists(a.file):
        print(f"{a.file} not found"); return
    raw = json.load(open(a.file))

    # run_robustness.py writes {"summary": {...}, "per_seed": [ {...}, ... ]},
    # window2d.py / ecl_sweep.py write a flat {key: run} mapping.  Accept both.
    summary = None
    if isinstance(raw, dict) and "per_seed" in raw:
        summary = raw.get("summary", {})
        raw = {str(r.get("seed", i)): r for i, r in enumerate(raw["per_seed"])}
        bits = [f"{k}={summary[k]}" for k in
                ("mode", "shared_frac", "balanced", "t_build", "t_gov")
                if k in summary]
        if bits:
            print("robustness file:", ", ".join(bits))
    if raw and not any(d.get("w") for d in raw.values()) \
            and any(d.get("A10") is not None for d in raw.values()):
        print("NOTE: no per-arm weights in this file; using the arm-acquisition")
        print("      probes A10/A01 (hub silenced, non-zero response) as the")
        print("      per-arm test of whether each arm cleared the cliff.\n")
    elif raw and not any(d.get("w_mean") is not None for d in raw.values()):
        print("WARNING: this file records neither arm weights nor arm-acquisition")
        print("         probes, so whether the arms were taught cannot be checked")
        print("         here.  Runs are classified on xor_index alone; treat the")
        print("         'taught' part of the verdict as inherited, not measured.\n")
    elif raw and not any("w" in d for d in raw.values()):
        print("NOTE: this file stores only the mean arm weight, not the two arms,")
        print("      so the asymmetric class cannot be detected here.  Re-run with")
        print("      the current run_robustness.py to record both.\n")

    # window2d.json values lack w_mean but carry w as a list; both are handled.
    runs = []
    for k, d in raw.items():
        new, leaky = rescore(d)
        # window2d.json stores w as a list and omits xor_index; derive both so
        # that the per-group medians are populated for either file format.
        g = lambda k, dv=None: (dv if d.get(k) is None else d[k])
        w_mean = g("w_mean")
        if w_mean is None and d.get("w"):
            w_mean = float(np.mean(d["w"]))
        xi = g("xor_index")
        if xi is None:
            smin = min(g("R10", 0.0), g("R01", 0.0))
            xi = (smin - g("R11", 0.0)) / smin if smin > 0 else 0.0
        runs.append(dict(key=k, old=d.get("cls", "?"), new=new, leaky=leaky,
                         grp=group_key(k, d), w_mean=w_mean, xor_index=xi,
                         **{f: g(f, 0.0) for f in ("R10", "R01", "R11", "R00")}))

    changed = [r for r in runs if r["old"] != r["new"]]
    print(f"{a.file}: {len(runs)} runs, {len(changed)} relabelled\n")

    if changed:
        tally = Counter((r["old"], r["new"]) for r in changed)
        print("relabelling summary (old -> new):")
        for (o, n), c in tally.most_common():
            print(f"  {o:>10s} -> {n:<10s}  {c:4d}")
        mis = [r for r in changed
               if r["old"] == "collapse" and r["new"] in ("XOR", "OR")]
        if mis:
            xis = [r["xor_index"] for r in mis if r["xor_index"] is not None]
            print(f"\n  {len(mis)} runs called 'collapse' were taught circuits.")
            if xis:
                print(f"  their xor_index: median {np.median(xis):.2f}, "
                      f"max {max(xis):.2f}")
            print(f"  (all had min(R10,R01) < {OLD_THETA_R:.0f} Hz, "
                  f"which is what the old gate keyed on)")

    grp = defaultdict(list)
    for r in runs:
        grp[r["grp"]].append(r)

    print("\nper-group, new labels (rate shown as descriptor only):")
    hdr = (f"{'group':>22} | {'n':>3} | {'XOR':>5} {'OR':>5} {'n_taught':>8} "
           f"{'asym':>5} {'silent':>6} | {'w med':>6} | {'xi med':>6} | {'R10 med':>8}")
    print(hdr); print("-" * len(hdr))
    for g in sorted(grp, key=lambda t: (str(t[0]), -(t[1] if t[1] is not None else 0))):
        rs = grp[g]
        n = len(rs)
        f = lambda lab: sum(r["new"] == lab for r in rs) / n
        wv = [r["w_mean"] for r in rs if r["w_mean"] is not None]
        xv = [r["xor_index"] for r in rs if r["xor_index"] is not None]
        rv = [r["R10"] for r in rs if r["R10"] is not None]
        print(f"{str(g):>22} | {n:3d} | {f('XOR'):5.2f} {f('OR'):5.2f} "
              f"{f('not_taught'):8.2f} {f('asymmetric'):5.2f} {f('silent'):6.2f} | "
              f"{np.median(wv) if wv else float('nan'):6.2f} | "
              f"{np.median(xv) if xv else float('nan'):6.2f} | "
              f"{np.median(rv) if rv else float('nan'):8.1f}")

    nleak = sum(r["leaky"] for r in runs)
    if nleak:
        print(f"\n{nleak} runs flagged leaky (R00 > {LEAK:.0f} Hz)")

    if a.csv:
        import csv
        with open(a.csv, "w", newline="") as fh:
            wtr = csv.DictWriter(fh, fieldnames=list(runs[0].keys()))
            wtr.writeheader()
            for r in runs:
                wtr.writerow(r)
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    main()
