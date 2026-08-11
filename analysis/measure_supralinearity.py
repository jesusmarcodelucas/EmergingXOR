#!/usr/bin/env python3
"""
measure_supralinearity.py -- the missing measurement for EmergingXOR.

Measures, in the circuit of phase_sim3.py, the central quantitative prediction
of the manuscript, which no existing script tested:

    (weak, Sec. 2)   g_inh(2) > 2 g_inh(1)
    (Eq. 2, Sec. 3)  g_inh(2) - 2 g_inh(1) > sigma_eff,   sigma_eff = gL + ge

Nothing is trained: arms are clamped at a fixed weight and only the hub->lateral
shunt is measured, so the result is a property of the architecture and its
operating point, not of one learning run.

It also runs the over-growth control that separates two mechanisms which the
manuscript previously conflated:
  - does 'builder only -> OR' come from over-grown arms, or from ECl never
    switching?  Probing the over-grown arms at MATURE ECl answers it.

    python3 measure_supralinearity.py            # both blocks
    python3 measure_supralinearity.py --quick    # shorter probes
"""
import sys as _sys, pathlib as _pathlib
_ROOT = _pathlib.Path(__file__).resolve().parents[1]
for _d in ("simulations", "analysis"):
    _p = str(_ROOT / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import argparse, time
import numpy as np
from brian2 import (NeuronGroup, PoissonGroup, Synapses, SpikeMonitor,
                    StateMonitor, Network, prefs, defaultclock, seed,
                    BrianLogger, ms, mV, nS, pA, Hz, second)
import two_phase_model as P

prefs.codegen.target = "numpy"
BrianLogger.suppress_name("resolution_conflict", "method_choice")
defaultclock.dt = 0.1 * ms


def probe(pattern, w_arm=0.45, w_hl=None, ecl=None, shared=1.0,
          seed_=11, T=0.6 * second):
    """Mean gi/ge in lateral 0 plus hub, lateral and read-out rates."""
    w_hl = P.W_HL_PER_SYN if w_hl is None else w_hl
    ecl = P.ECL_MATURE if ecl is None else ecl
    seed(seed_)
    ns = dict(Cm=P.Cm, gL=P.gL, EL=P.EL, VT=P.VT, Vr=P.Vr, Ee=P.Ee,
              tau_e=P.tau_e, tau_i=P.tau_i, tau_pre=P.tau_pre, tau_Ca=P.tau_Ca,
              kappa=P.kappa, tau_theta=P.tau_theta, W_SCALE=P.W_SCALE, eta=0.0,
              tau_w=P.tau_w, b_adapt=P.b_adapt, WMIN=P.WMIN, WMAX=P.WMAX)
    Inp = PoissonGroup(2, rates=np.asarray(pattern) * P.IN_RATE)
    Lat = NeuronGroup(2, P.lat_eqs, threshold="v>VT",
                      reset="v=Vr; Ca+=1; wad+=b_adapt", refractory=2 * ms,
                      method="euler", namespace=ns)
    Lat.v = P.EL; Lat.wad = 0 * pA; Lat.theta = 0.35; Lat.ECl = ecl
    Hub = NeuronGroup(P.NHUB, P.hub_eqs, threshold="v>VT_h", reset="v=Vr",
                      refractory=2 * ms, method="euler", namespace=ns)
    Hub.v = P.EL
    Hub.VT_h = P.VT + np.linspace(-1.0, 1.0, P.NHUB) * P.VT_H_SPREAD
    Read = NeuronGroup(1, P.simple_eqs, threshold="v>VT", reset="v=Vr",
                       refractory=2 * ms, method="euler", namespace=ns)
    Read.v = P.EL
    Arm = Synapses(Inp, Lat, model="""w : 1
                                      dx/dt = -x/tau_pre : 1 (clock-driven)""",
                   on_pre="ge_post += w*W_SCALE; x += 1", method="euler",
                   namespace=ns)
    Arm.connect(j="i"); Arm.w = w_arm
    S_ih = Synapses(Inp, Hub, on_pre="ge_post += w_h",
                    namespace=dict(w_h=P.W_HUB * shared)); S_ih.connect(True)
    S_hl = Synapses(Hub, Lat, on_pre="gi_post += w_sh",
                    namespace=dict(w_sh=w_hl)); S_hl.connect(True)
    S_lr = Synapses(Lat, Read, on_pre="ge_post += W_READ",
                    namespace=dict(W_READ=P.W_READ)); S_lr.connect(True)
    gm = StateMonitor(Lat, ["gi", "ge"], record=[0], dt=1 * ms)
    hs, rs = SpikeMonitor(Hub), SpikeMonitor(Read)
    Network([Inp, Lat, Hub, Read, Arm, S_ih, S_hl, S_lr, gm, hs, rs]).run(T)
    Ts = float(T / second); sk = 200
    return dict(gi=float(np.mean(gm.gi[0][sk:] / nS)),
                ge=float(np.mean(gm.ge[0][sk:] / nS)),
                hub=hs.num_spikes / (P.NHUB * Ts), read=rs.count[0] / Ts)


def row(label, T, **kw):
    o = {k: probe(p, T=T, **kw) for k, p in
         [("10", (1, 0)), ("01", (0, 1)), ("11", (1, 1))]}
    g1 = 0.5 * (o["10"]["gi"] + o["01"]["gi"]); g2 = o["11"]["gi"]
    ratio = g2 / g1 if g1 > 1e-9 else float("nan")
    excess = g2 - 2 * g1; sigma = float(P.gL / nS) + o["11"]["ge"]
    smin = min(o["10"]["read"], o["01"]["read"])
    xi = (smin - o["11"]["read"]) / smin if smin > 0 else float("nan")
    cls = "collapse" if smin < 25 else ("OR" if o["11"]["read"] > 0.5 * smin else "XOR")
    print(f"{label:>14} | hub {o['10']['hub']:6.1f}->{o['11']['hub']:6.1f}Hz "
          f"| g_inh {g1:6.1f}->{g2:6.1f}nS ({ratio:5.2f}x) "
          f"{'SUPRA' if ratio > 2 else '  sub'} | excess {excess:6.1f} vs "
          f"sigma {sigma:5.1f} {'OK' if excess > sigma else 'no'} "
          f"| R {o['10']['read']:6.1f}/{o['11']['read']:6.1f} xi={xi:5.2f} {cls}",
          flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    T = (0.4 if a.quick else 0.6) * second
    t0 = time.time()
    print(f"phase_sim3: IN_RATE={P.IN_RATE/Hz:.0f}Hz W_HUB={P.W_HUB/nS:.1f}nS "
          f"NHUB={P.NHUB} W_HL={P.W_HL_PER_SYN/nS:.0f}nS "
          f"ECl_mature={P.ECL_MATURE/mV:.0f}mV\n")
    print("== supralinearity of hub recruitment (arms clamped) ==")
    for w in [0.10, 0.45, 0.77]:
        row(f"arm w={w:.2f}", T, w_arm=w)
    print("\n== over-growth control: probed at MATURE ECl, hub active ==")
    print("   (does over-growth alone give OR, or is 'builder only -> OR' an")
    print("    artefact of probing with ECl still depolarising?)")
    for w in [0.45, 0.77, 1.20, 1.80, 2.60]:
        row(f"arm w={w:.2f}", T, w_arm=w)
    print(f"\n[{time.time()-t0:.0f}s]")
