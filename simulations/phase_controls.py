#!/usr/bin/env python3
"""
run_phase_experiments.py  --  run the hardened two-phase model (phase_sim2.py) and
write the data files that make_phase_figures.py turns into figures:

    phase_factorial.json   : the 5 factorial conditions + arm-acquisition verdicts
    phase_traces2.npz      : clean recorded traces of the FULL sequence

Reduced durations (t_build=18, t_ramp=4, t_gov=12) are used so the whole thing
runs in a few minutes; they sit inside the XOR window.  Increase to
t_build=25, t_gov=25 to reproduce the "canonical" headline numbers.
Runtime is dominated by Brian2; expect a few minutes with the numpy target.
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import json, numpy as np
import two_phase_model as P


def main():
    TB, TR, TG, SEED = 18.0, 4.0, 12.0, 11

    rows = []
    for mode, label in [("governor_only", "governor only"),
                        ("no_hub", "no hub"),
                        ("builder_only", "builder only"),
                        ("full", "FULL sequence")]:
        r = P.run_phases(mode=mode, t_build=TB, t_ramp=TR, t_gov=TG, seed_=SEED)
        rows.append(dict(label=label,
                         vals=[round(r["R10"]), round(r["R01"]), round(r["R11"]), round(r["R00"])],
                         cls=P.classify(r), arms=f"[{r['w'][0]:.2f}, {r['w'][1]:.2f}]",
                         acq=P.arm_verdict(r)))
        print(f"{label:<15} {rows[-1]['vals']} -> {rows[-1]['cls']:8s} | arms {rows[-1]['acq']}", flush=True)

    # eta=0 control (full transition, plasticity frozen)
    r0 = P.run_phases(mode="full", eta=0.0, t_build=TB, t_ramp=TR, t_gov=TG, seed_=SEED)
    eta0 = dict(label="no plasticity ($\\eta{=}0$)",
                vals=[round(r0["R10"]), round(r0["R01"]), round(r0["R11"]), round(r0["R00"])],
                cls=P.classify(r0), arms=f"[{r0['w'][0]:.2f}, {r0['w'][1]:.2f}]", acq=P.arm_verdict(r0))
    print(f"{'eta=0':<15} {eta0['vals']} -> {eta0['cls']:8s} | arms {eta0['acq']}", flush=True)

    # order for the figure: governor, no hub, builder, eta0, full
    ordered = [rows[0], rows[1], rows[2], eta0, rows[3]]
    json.dump(ordered, open("phase_factorial.json", "w"), indent=1)
    print("wrote phase_factorial.json")

    # clean recorded traces of the full sequence
    rr = P.run_phases(mode="full", t_build=TB, t_ramp=TR, t_gov=TG, seed_=SEED, record=True)
    np.savez("phase_traces2.npz", t_w=rr["t_w"], W=rr["W"], t_l=rr["t_l"],
             TH=rr["TH"], C=rr["C"], ECL=rr["ECL"])
    print("wrote phase_traces2.npz  (arms %.2f->%.2f)" % (float(rr["W"].min()), float(rr["W"].max())))


if __name__ == "__main__":
    main()
