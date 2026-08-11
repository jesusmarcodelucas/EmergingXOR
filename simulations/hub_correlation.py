#!/usr/bin/env python3
"""
run_robustness.py -- second-round robustness battery for phase_sim3.py.

Runs MANY seeds, classifies EACH run individually (never the average), and
reports the distribution the reviewer asked for:
  - XOR fraction, both-grown fraction, asymmetric fraction
  - median / IQR of R10, R01, R11, xor_index, symmetry

It also sweeps `hub_shared_frac` (how correlated the hub afferents are), because
a short exploratory run showed the XOR result is NOT robust to partially
independent hub inputs: at shared_frac=0.5 the circuit computes OR, not XOR.
This battery is meant to locate the boundary and quantify it with statistics.

USAGE (parallelise across your DGX cores):

    # baseline, 100 seeds, fully-shared hub (the published condition)
    python3 run_robustness.py --seeds 100 200 --shared 1.0 --workers 32 \
        --t_build 18 --t_gov 12 --out robust_shared1.00.json

    # hub-independence sweep, 40 seeds each
    for sf in 1.0 0.95 0.9 0.8 0.7 0.6 0.5; do
      python3 run_robustness.py --seeds 100 140 --shared $sf --workers 32 \
          --t_build 18 --t_gov 12 --out robust_shared${sf}.json
    done

    # fully-random (unbalanced) pattern presentation, as a supplementary control
    python3 run_robustness.py --seeds 100 140 --shared 1.0 --balanced 0 --workers 32 \
        --out robust_unbalanced.json

Use --t_build 25 --t_gov 25 for the canonical durations (slower).
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import argparse, json, os
from multiprocessing import Pool
import numpy as np


def _one(job):
    # import inside the worker so each process initialises Brian2 cleanly
    import two_phase_model as P
    mode, seed_, tb, tr, tg, bal, sf = job
    r = P.run_phases(mode=mode, seed_=seed_, t_build=tb, t_ramp=tr, t_gov=tg,
                     balanced=bool(bal), hub_shared_frac=sf)
    keep = ("seed", "R10", "R01", "R11", "R00", "A10", "A01", "xor_index", "symmetry")
    out = {k: float(r[k]) for k in keep}
    out["cls"] = P.classify(r)
    out["arm"] = P.arm_verdict(r)
    out["w_mean"] = float(np.mean(r["w"]))
    out["w"] = [float(x) for x in r["w"]]   # both arms: needed to detect asymmetry
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="full")
    ap.add_argument("--seeds", type=int, nargs=2, default=[100, 200])
    ap.add_argument("--shared", type=float, default=1.0)
    ap.add_argument("--balanced", type=int, default=1)
    ap.add_argument("--t_build", type=float, default=18.0)
    ap.add_argument("--t_ramp", type=float, default=4.0)
    ap.add_argument("--t_gov", type=float, default=12.0)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.out is None:
        a.out = f"robust_{a.mode}_shared{a.shared:.2f}_bal{a.balanced}.json"
    seeds = list(range(a.seeds[0], a.seeds[1]))
    jobs = [(a.mode, s, a.t_build, a.t_ramp, a.t_gov, a.balanced, a.shared) for s in seeds]
    print(f"{len(jobs)} runs | mode={a.mode} shared_frac={a.shared} balanced={a.balanced} "
          f"t_build={a.t_build} | {a.workers} workers", flush=True)

    runs = []
    with Pool(a.workers) as pool:
        for k, out in enumerate(pool.imap_unordered(_one, jobs), 1):
            runs.append(out)
            print(f"[{k:3d}/{len(jobs)}] seed {out['seed']:.0f}: {out['cls']:8s} | "
                  f"arm {out['arm']:12s} | {out['R10']:.0f}/{out['R01']:.0f}/{out['R11']:.0f} | "
                  f"xor_idx={out['xor_index']:.2f} sym={out['symmetry']:.2f}", flush=True)

    def iqr(key):
        v = np.array([r[key] for r in runs])
        return [float(np.percentile(v, 50)), float(np.percentile(v, 25)), float(np.percentile(v, 75))]

    S = dict(
        n=len(runs), mode=a.mode, shared_frac=a.shared, balanced=a.balanced,
        t_build=a.t_build, t_gov=a.t_gov,
        xor_fraction=float(np.mean([r["cls"] == "XOR" for r in runs])),
        or_fraction=float(np.mean([r["cls"] == "OR" for r in runs])),
        collapse_fraction=float(np.mean([r["cls"] == "collapse" for r in runs])),
        both_grown_fraction=float(np.mean([r["arm"] == "both grown" for r in runs])),
        asymmetric_fraction=float(np.mean([r["arm"] == "asymmetric" for r in runs])),
        R10=iqr("R10"), R01=iqr("R01"), R11=iqr("R11"),
        xor_index=iqr("xor_index"), symmetry=iqr("symmetry"),
        w_mean_median=float(np.median([r["w_mean"] for r in runs])),
    )
    print("\n==== SUMMARY ====")
    for k in ("n", "xor_fraction", "or_fraction", "collapse_fraction",
              "both_grown_fraction", "asymmetric_fraction"):
        print(f"  {k:20s}: {S[k]}")
    for k in ("R10", "R01", "R11", "xor_index", "symmetry"):
        print(f"  {k:20s}: median {S[k][0]:.2f}  IQR [{S[k][1]:.2f}, {S[k][2]:.2f}]")

    os.makedirs("data", exist_ok=True)
    out = a.out or f"data/robust_shared{a.shared:.2f}_tb{int(a.t_build)}_{a.mode}_bal{a.balanced}.json"
    json.dump({"summary": S, "per_seed": runs}, open(out, "w"), indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
