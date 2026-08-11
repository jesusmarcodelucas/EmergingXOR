#!/usr/bin/env python3
"""
reproduce_all.py -- regenerate every figure in the preprint.

Stages run in dependency order: simulations write caches into data/, figure
scripts read those caches. Each stage is skipped if its output already exists,
so an interrupted run resumes rather than restarting.

    python reproduce_all.py --quick     # schematics + cached figures (~10 min)
    python reproduce_all.py             # everything, several hours
    python reproduce_all.py --check     # verify consistency, run nothing
    python reproduce_all.py --list      # show the plan and exit

--quick is the one to run first. If it fails, the environment is wrong and
there is no point launching the long sweeps.
"""
import argparse
import json
import os
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
OUT = ROOT / "figures" / "output"

# (label, argv, produced-file, long?)   paths relative to ROOT
STAGES = [
    # Fig. 1 is TikZ, inlined in the manuscript: nothing to generate here.
    # See figures/tikz/fig1_mature_motif.tex and its standalone preview.
    ("Fig. 2  developmental sequence (schematic)",
     ["figures/fig2_development.py"], "motif_dev_en.png", False),
    ("Fig. 9  builder-governor (schematic)",
     ["figures/fig9_builder_governor.py"], "preprint_fig_en.png", False),
    ("Fig. 3  GABA switch network",
     ["simulations/gaba_switch.py"], "demo_brian2_en.png", False),
    ("Figs. 5-6  data: growth traces and lesion controls",
     ["simulations/phase_controls.py"], "data/fig5_6/phase_factorial.json", False),
    ("Figs. 5-6  figures",
     ["figures/fig5_fig6_growth_controls.py"], "phase_trace_en.png", False),
    ("Fig. 4  data: chloride set-point sweep",
     ["simulations/ecl_sweep.py", "--seeds", "20"], "data/fig4/ecl_sweep.json", True),
    ("Fig. 4  figure",
     ["figures/fig4_ecl_setpoints.py"], "ecl_sweep_en.png", False),
    ("Fig. 7  data: developmental window",
     ["simulations/developmental_window.py", "--seeds", "20"],
     "data/fig7/window2d.json", True),
    ("Fig. 7  figure",
     ["simulations/developmental_window.py", "--plot"], "window2d_en.png", False),
    ("Fig. 8  figure (needs the hub-correlation sweep, see README)",
     ["figures/fig8_hub_correlation.py"], "robustness_sweep_en.png", False),
]


def produced(rel):
    """Figures land in the cwd of their script; caches have explicit paths."""
    p = pathlib.Path(rel)
    if p.suffix == ".png":
        for cand in (OUT / p.name, ROOT / p.name,
                     ROOT / "figures" / p.name, ROOT / "simulations" / p.name):
            if cand.exists():
                return cand
        return None
    return (ROOT / p) if (ROOT / p).exists() else None


def check():
    """Assert that the two copies of the classification rule agree."""
    sys.path.insert(0, str(ROOT / "analysis"))
    sys.path.insert(0, str(ROOT / "simulations"))
    import rescore_runs
    import developmental_window as dw

    cache = ROOT / "data/fig7/window2d.json"
    if not cache.exists():
        print("[check] data/fig7/window2d.json absent; nothing to cross-check.")
        return True
    runs = json.load(open(cache))
    bad = []
    for k, d in runs.items():
        a = rescore_runs.rescore(d)[0]
        b = dw.label(d)
        if a != b:
            bad.append((k, a, b))
    if bad:
        print(f"[check] MISMATCH between rescore_runs.rescore() and "
              f"developmental_window.label() on {len(bad)}/{len(runs)} runs:")
        for k, a, b in bad[:5]:
            print(f"          {k}: rescore={a!r}  window={b!r}")
        return False
    print(f"[check] classification agrees on all {len(runs)} cached runs.")
    return True


def run(stage, dry):
    label, argv, out, _ = stage
    have = produced(out)
    if have:
        print(f"[skip] {label}\n       exists: {have.relative_to(ROOT)}")
        return True
    if dry:
        print(f"[todo] {label}\n       -> {out}")
        return True
    print(f"[run ] {label}", flush=True)
    t0 = time.time()
    script = ROOT / argv[0]
    cwd = script.parent
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT / "simulations"), str(ROOT / "analysis"),
         env.get("PYTHONPATH", "")])
    r = subprocess.run([sys.executable, str(script.name)] + argv[1:],
                       cwd=cwd, env=env)
    dt = time.time() - t0
    if r.returncode != 0:
        print(f"[FAIL] {label}  (exit {r.returncode}, {dt/60:.1f} min)")
        return False
    print(f"[ok  ] {label}  ({dt/60:.1f} min)")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip the long sweeps (Figs. 4 and 7 data)")
    ap.add_argument("--list", action="store_true", help="show the plan, run nothing")
    ap.add_argument("--check", action="store_true",
                    help="verify classification consistency, run nothing")
    a = ap.parse_args()

    if a.check:
        sys.exit(0 if check() else 1)

    stages = [s for s in STAGES if not (a.quick and s[3])]
    if a.quick:
        print("quick mode: long sweeps skipped, cached data used where present\n")

    ok = True
    for s in stages:
        if not run(s, dry=a.list):
            ok = False
            break

    if not a.list and ok:
        print()
        ok = check()
    print("\n" + ("done" if ok else "stopped on error"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
