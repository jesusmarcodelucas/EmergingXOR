#!/usr/bin/env python3
"""The honest Fig. 5 for the two-phase model.

Sweeps the TIMING of the switch against the STRENGTH of the mature shunt, and
asks how wide the XOR window really is. Resumable; parallel.

    python3 window2d.py --seeds 3 --workers 8
"""
import argparse, json, os
from collections import Counter
from multiprocessing import Pool

T_BUILD = [8.0, 15.0, 20.0, 25.0, 32.0, 40.0]      # when the switch happens (s)
W_HL    = [10.0, 16.0, 22.0, 28.0, 34.0]           # mature shunt, nS per synapse
CACHE   = "window2d.json"


def one(job):
    from brian2 import nS
    import phase_sim as P
    tb, whl, sd = job
    r = P.run_phases(mode="full", t_build=tb, w_hl=whl * nS, seed_=sd)
    return (f"{tb}|{whl}|{sd}",
            dict(cls=P.classify(r), R10=r["R10"], R01=r["R01"],
                 R11=r["R11"], R00=r["R00"], w=[float(x) for x in r["w"]]))


def save(res):
    with open(CACHE + ".tmp", "w") as fh:
        json.dump(res, fh, indent=1)
    os.replace(CACHE + ".tmp", CACHE)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    seeds = [11, 23, 37, 51, 67][:a.seeds]

    res = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    jobs = [(tb, whl, sd) for tb in T_BUILD for whl in W_HL for sd in seeds
            if f"{tb}|{whl}|{sd}" not in res]
    total = len(T_BUILD) * len(W_HL) * len(seeds)
    print(f"{total} runs, {len(jobs)} to do, {a.workers} workers", flush=True)

    with Pool(a.workers) as pool:
        for k, v in pool.imap_unordered(one, jobs):
            res[k] = v
            save(res)
            print(f"[{len(res):3d}/{total}] {k} -> {v['cls']:8s} "
                  f"{v['R10']:5.0f}/{v['R01']:5.0f}/{v['R11']:5.0f}", flush=True)

    print("\n" + "=" * 62)
    print(f"{'t_build':>8} |" + "".join(f"{w:>10.0f}" for w in W_HL))
    print("-" * 62)
    for tb in T_BUILD:
        row = f"{tb:8.0f} |"
        for whl in W_HL:
            outs = [res[f"{tb}|{whl}|{sd}"]["cls"] for sd in seeds
                    if f"{tb}|{whl}|{sd}" in res]
            if outs:
                dom, n = Counter(outs).most_common(1)[0]
                row += f"{dom[:3].upper() + f' {n}/{len(outs)}':>10}"
        print(row)
    print("\nXOR cells:", sum(v["cls"] == "XOR" for v in res.values()),
          "/", len(res))
