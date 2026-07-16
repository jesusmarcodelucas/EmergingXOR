#!/usr/bin/env python3
"""Is there a developmental WINDOW for the switch?

Too early : the arms never clear the firing cliff              -> collapse
Too late  : the builder over-grows them and no shunt can then
            suppress coincidence                               -> OR
In between:                                                    -> XOR
"""
from collections import Counter
import numpy as np
import phase_sim as P

T_BUILD = [5.0, 12.0, 18.0, 25.0, 32.0, 40.0, 50.0]
SEEDS = [11, 23, 37]

print("switch-timing sweep (arms start BELOW the cliff, w0 = 0.10)\n")
print(f"{'t_build':>8} | {'arms':>6} | {'(1,0)':>6} {'(0,1)':>6} {'(1,1)':>6} | outcome (n=3)")
print("-" * 64)
for tb in T_BUILD:
    outs, arms, a, b, c = [], [], [], [], []
    for sd in SEEDS:
        r = P.run_phases(mode="full", t_build=tb, seed_=sd)
        outs.append(P.classify(r)); arms.append(r["w"].mean())
        a.append(r["R10"]); b.append(r["R01"]); c.append(r["R11"])
    dom, n = Counter(outs).most_common(1)[0]
    print(f"{tb:8.0f} | {np.mean(arms):6.2f} | {np.mean(a):6.0f} "
          f"{np.mean(b):6.0f} {np.mean(c):6.0f} | {dom} {n}/3", flush=True)
