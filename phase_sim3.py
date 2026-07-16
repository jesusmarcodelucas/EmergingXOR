#!/usr/bin/env python3
"""
phase_sim3.py -- two-phase model, second-round-hardened.

Adds, on top of phase_sim2.py (which fixed the state-reset / recording / eta=0 /
arm-acquisition / matched-RNG issues):

  * MULTI-SEED harness that classifies EACH run individually and reports the
    distribution (XOR fraction, both-grown fraction, median/IQR of rates,
    arm-weight distribution, asymmetric failures) -- never classifies averaged rates.
  * BALANCED shuffled pattern block: every 4 presentations contain each of the
    four input patterns exactly once, order random (fully-random still available).
  * Total hub strength parameterised (W_HL_TOTAL) so changing NHUB does not
    silently change the total inhibitory conductance onto a lateral.
  * PARTIALLY-INDEPENDENT hub afferents: `hub_shared_frac` in [0,1]. At 1.0 all
    hub cells share the two Poisson trains (as before); below 1.0 each hub cell
    also receives its own independent pattern-driven Poisson source, so the
    population is no longer maximally correlated. Robustness knob.
  * Explicit probe signature probe(pattern, hub_to_lat=, ecl=) -- the
    arm-acquisition control is now unambiguous.
  * ECl recorded in its own monitor (record=[0]); no shared-variable warning.
  * Continuous indices: xor_index = (min_single - R11)/min_single ; symmetry.
  * The default 'full' condition is computed once, not twice.
"""
import numpy as np
from brian2 import *

prefs.codegen.target = "numpy"          # set "cython" for speed if available
BrianLogger.suppress_name("resolution_conflict", "method_choice")
defaultclock.dt = 0.1 * ms

Cm, gL = 200 * pF, 10 * nS
EL, VT, Vr = -65 * mV, -50 * mV, -60 * mV
Ee = 0 * mV
tau_e, tau_i = 5 * ms, 10 * ms
tau_w, b_adapt = 250 * ms, 90 * pA
ECL_IMMATURE, ECL_MATURE = -40 * mV, -65 * mV

tau_pre, tau_Ca = 15 * ms, 60 * ms
kappa, tau_theta = 9.0, 800 * ms
ETA_STEP = 3e-5
WMIN, WMAX = 0.05, 2.6

IN_RATE = 200 * Hz
NHUB = 20
W_SCALE = 14.2 * nS
W_HUB = 2.2 * nS
W_HL_TOTAL = 440.0 * nS                  # total hub->lateral conductance ...
W_HL_PER_SYN = W_HL_TOTAL / NHUB         # ... normalised by NHUB (= 22 nS at NHUB=20)
W_READ = 60.0 * nS
W0_BELOW_CLIFF = 0.10
VT_H_SPREAD = 1.5 * mV

lat_eqs = """
dv/dt   = (gL*(EL-v) + ge*(Ee-v) + gi*(ECl-v) - wad)/Cm : volt (unless refractory)
dge/dt  = -ge/tau_e : siemens
dgi/dt  = -gi/tau_i : siemens
dwad/dt = -wad/tau_w : amp
dCa/dt  = -Ca/tau_Ca : 1
c       = Ca/kappa : 1
dtheta/dt = (c*c - theta)/tau_theta : 1
ECl     : volt (shared)
"""
hub_eqs = """
dv/dt  = (gL*(EL-v) + ge*(Ee-v))/Cm : volt (unless refractory)
dge/dt = -ge/tau_e : siemens
VT_h   : volt
"""
simple_eqs = """
dv/dt  = (gL*(EL-v) + ge*(Ee-v))/Cm : volt (unless refractory)
dge/dt = -ge/tau_e : siemens
"""
PATTERNS = [(1, 0), (0, 1), (1, 1), (0, 0)]
VALID = {"full", "governor_only", "builder_only", "no_hub"}


def run_phases(mode="full", w0=W0_BELOW_CLIFF, w_hl=W_HL_PER_SYN, eta=ETA_STEP,
               t_build=25.0, t_ramp=5.0, t_gov=25.0, seed_=11, record=False,
               balanced=True, hub_shared_frac=1.0):
    if mode not in VALID:
        raise ValueError(f"Unknown mode {mode!r}; expected {sorted(VALID)}")
    if min(t_build, t_ramp, t_gov) < 0:
        raise ValueError("Phase durations must be non-negative")
    if not WMIN <= w0 <= WMAX:
        raise ValueError(f"w0 must be in [{WMIN},{WMAX}]")
    if not 0.0 <= hub_shared_frac <= 1.0:
        raise ValueError("hub_shared_frac must be in [0,1]")
    seed(seed_)
    rng = np.random.default_rng(seed_ + 1)
    net = Network()
    ns = dict(Cm=Cm, gL=gL, EL=EL, VT=VT, Vr=Vr, Ee=Ee, tau_e=tau_e, tau_i=tau_i,
              tau_pre=tau_pre, tau_Ca=tau_Ca, kappa=kappa, tau_theta=tau_theta,
              W_SCALE=W_SCALE, eta=eta, tau_w=tau_w, b_adapt=b_adapt,
              WMIN=WMIN, WMAX=WMAX)

    Inp = PoissonGroup(2, rates=[0, 0] * Hz)
    Lat = NeuronGroup(2, lat_eqs, threshold="v>VT",
                      reset="v=Vr; Ca+=1; wad+=b_adapt", refractory=2 * ms,
                      method="euler", namespace=ns)
    Lat.v = EL; Lat.wad = 0 * pA; Lat.theta = 0.35
    Lat.ECl = ECL_MATURE if mode == "governor_only" else ECL_IMMATURE

    Hub = NeuronGroup(NHUB, hub_eqs, threshold="v>VT_h", reset="v=Vr",
                      refractory=2 * ms, method="euler", namespace=ns)
    Hub.v = EL + 5 * mV * np.random.rand(NHUB)
    Hub.VT_h = VT + np.linspace(-1.0, 1.0, NHUB) * VT_H_SPREAD
    Read = NeuronGroup(1, simple_eqs, threshold="v>VT", reset="v=Vr",
                       refractory=2 * ms, method="euler", namespace=ns)
    Read.v = EL

    Arm = Synapses(Inp, Lat,
                   model="""w : 1
                            dx/dt = -x/tau_pre : 1 (clock-driven)""",
                   on_pre="ge_post += w*W_SCALE; x += 1", method="euler", namespace=ns)
    Arm.connect(j="i"); Arm.w = w0
    plast = Arm.run_regularly(
        "w = clip(w + eta*x*c_post*(c_post - theta_post), WMIN, WMAX)",
        dt=1 * ms, when="end")

    # hub afferents: shared component (all hubs see the two Inp trains) ...
    S_ih = Synapses(Inp, Hub, on_pre="ge_post += w_h",
                    namespace=dict(w_h=W_HUB * hub_shared_frac)); S_ih.connect(True)
    IndpInp = None
    if hub_shared_frac < 1.0:                 # ... plus per-hub independent trains
        IndpInp = PoissonGroup(2 * NHUB, rates=[0] * (2 * NHUB) * Hz)
        S_indp = Synapses(IndpInp, Hub, on_pre="ge_post += w_h",
                          namespace=dict(w_h=W_HUB * (1.0 - hub_shared_frac)))
        S_indp.connect(condition="int(i/2) == j")

    hl = 0 * nS if mode == "no_hub" else w_hl
    S_hl = Synapses(Hub, Lat, on_pre="gi_post += w_sh",
                    namespace=dict(w_sh=hl)); S_hl.connect(True)
    S_lr = Synapses(Lat, Read, on_pre="ge_post += W_READ",
                    namespace=dict(W_READ=W_READ)); S_lr.connect(True)
    net.add(Inp, Lat, Hub, Read, Arm, S_ih, S_hl, S_lr)
    if IndpInp is not None:
        net.add(IndpInp, S_indp)

    queue = []

    @network_operation(dt=200 * ms)
    def present():
        if balanced:
            if not queue:
                queue.extend(rng.permutation(len(PATTERNS)).tolist())
            pat = PATTERNS[queue.pop()]
        else:
            pat = PATTERNS[rng.integers(len(PATTERNS))]
        r = np.asarray(pat) * IN_RATE
        Inp.rates = r
        if IndpInp is not None:
            IndpInp.rates = np.tile(r, NHUB)
    net.add(present)

    mons = {}
    if record:
        mons["l"] = StateMonitor(Lat, ["theta", "c"], record=True, dt=25 * ms)
        mons["ecl"] = StateMonitor(Lat, "ECl", record=[0], dt=25 * ms)
        mons["w"] = StateMonitor(Arm, "w", record=True, dt=50 * ms)
        for m in mons.values():
            net.add(m)

    if mode != "governor_only":
        net.run(t_build * second)
    if mode not in ("governor_only", "builder_only") and t_ramp > 0:
        rs = float(net.t / second)
        @network_operation(dt=50 * ms, when="start")
        def ramp():
            f = np.clip((float(net.t / second) - rs) / t_ramp, 0.0, 1.0)
            Lat.ECl = ECL_IMMATURE + f * (ECL_MATURE - ECL_IMMATURE)
        net.add(ramp); net.run(t_ramp * second); net.remove(ramp)
        Lat.ECl = ECL_MATURE
    if mode == "builder_only":
        net.run((t_ramp + t_gov) * second)
    elif mode == "governor_only":
        Lat.ECl = ECL_MATURE; net.run((t_build + t_ramp + t_gov) * second)
    else:
        Lat.ECl = ECL_MATURE; net.run(t_gov * second)

    w_trained = np.array(Arm.w[:]).copy()
    theta_trained = np.array(Lat.theta[:]).copy()

    recorded = None
    if record:
        recorded = dict(t_w=np.asarray(mons["w"].t / second).copy(),
                        W=np.asarray(mons["w"].w).copy(),
                        t_l=np.asarray(mons["l"].t / second).copy(),
                        TH=np.asarray(mons["l"].theta).copy(),
                        C=np.asarray(mons["l"].c).copy(),
                        ECL=np.asarray(mons["ecl"].ECl[0] / mV).copy())
        for m in mons.values():
            m.active = False

    plast.active = False
    net.remove(present)
    if IndpInp is not None:
        IndpInp.rates = [0] * (2 * NHUB) * Hz
    probe_ecl = ECL_IMMATURE if mode == "builder_only" else ECL_MATURE

    def hard_reset():
        Inp.rates = 0 * Hz
        Lat.v = EL; Lat.ge = 0 * nS; Lat.gi = 0 * nS; Lat.wad = 0 * pA; Lat.Ca = 0
        Hub.v = EL; Hub.ge = 0 * nS
        Read.v = EL; Read.ge = 0 * nS
        Arm.x = 0

    Lat.ECl = probe_ecl
    hard_reset(); net.run(50 * ms)
    Mr, Ml, Mh = SpikeMonitor(Read), SpikeMonitor(Lat), SpikeMonitor(Hub)
    net.add(Mr, Ml, Mh); net.store("probe")

    def probe(pattern, dur=1.0, hub_to_lat=True, ecl=None):
        net.restore("probe", restore_random_state=True)
        S_hl.active = hub_to_lat
        hard_reset()
        Lat.ECl = probe_ecl if ecl is None else ecl
        Inp.rates = np.asarray(pattern) * IN_RATE
        net.run(dur * second)
        return float(Mr.count[0]) / dur, np.asarray(Ml.count) / dur

    R10, _ = probe((1, 0)); R01, _ = probe((0, 1))
    R11, _ = probe((1, 1)); R00, _ = probe((0, 0), dur=0.5)
    # arm-acquisition: mature ECl, hub silenced -- isolates what the arms learned
    A10, _ = probe((1, 0), hub_to_lat=False, ecl=ECL_MATURE)
    A01, _ = probe((0, 1), hub_to_lat=False, ecl=ECL_MATURE)
    S_hl.active = True

    res = dict(mode=mode, w=w_trained, theta=theta_trained,
               R10=R10, R01=R01, R11=R11, R00=R00, A10=A10, A01=A01,
               hub_shared_frac=hub_shared_frac, seed=seed_)
    res.update(metrics(res))
    if recorded:
        res.update(recorded)
    return res


def classify(r, theta_R=25.0, alpha=0.5, theta_0=10.0):
    smin = min(r["R10"], r["R01"])
    if smin < theta_R:
        return "collapse"
    if r["R00"] > theta_0:
        return "leaky"
    if r["R11"] > alpha * smin:
        return "OR"
    return "XOR"


def arm_verdict(r, thr=25.0):
    a, b = r.get("A10", 0.0), r.get("A01", 0.0)
    if a >= thr and b >= thr:
        return "both grown"
    if a >= thr or b >= thr:
        return "asymmetric"
    return "neither grown"


def metrics(r):
    """Continuous, threshold-free indices."""
    s = min(r["R10"], r["R01"]); mx = max(r["R10"], r["R01"], 1e-9)
    return dict(xor_index=(s - r["R11"]) / max(s, 1e-9),   # 1 = perfect XOR, <=0 = OR/none
                symmetry=s / mx)                            # 1 = symmetric arms


def summarise(runs):
    """Aggregate a list of per-run result dicts -- classify EACH, then pool."""
    n = len(runs)
    cls = [classify(r) for r in runs]
    av = [arm_verdict(r) for r in runs]
    def iqr(key):
        v = np.array([r[key] for r in runs])
        return np.percentile(v, 50), np.percentile(v, 25), np.percentile(v, 75)
    out = dict(
        n=n,
        xor_fraction=np.mean([c == "XOR" for c in cls]),
        or_fraction=np.mean([c == "OR" for c in cls]),
        collapse_fraction=np.mean([c == "collapse" for c in cls]),
        both_grown_fraction=np.mean([a == "both grown" for a in av]),
        asymmetric_fraction=np.mean([a == "asymmetric" for a in av]),
        R10=iqr("R10"), R01=iqr("R01"), R11=iqr("R11"),
        xor_index=iqr("xor_index"), symmetry=iqr("symmetry"),
        w_median=float(np.median([np.mean(r["w"]) for r in runs])),
    )
    return out


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="full")
    ap.add_argument("--seeds", type=int, nargs=2, default=[100, 105],
                    help="seed range [start, stop)")
    ap.add_argument("--t_build", type=float, default=18.0)
    ap.add_argument("--t_gov", type=float, default=12.0)
    ap.add_argument("--t_ramp", type=float, default=4.0)
    ap.add_argument("--balanced", type=int, default=1)
    ap.add_argument("--hub_shared_frac", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    seeds = list(range(a.seeds[0], a.seeds[1]))
    runs = []
    for s in seeds:
        r = run_phases(mode=a.mode, seed_=s, t_build=a.t_build, t_ramp=a.t_ramp,
                       t_gov=a.t_gov, balanced=bool(a.balanced),
                       hub_shared_frac=a.hub_shared_frac)
        runs.append(r)
        print(f"seed {s}: {classify(r):8s} | arms {arm_verdict(r):12s} | "
              f"{r['R10']:.0f}/{r['R01']:.0f}/{r['R11']:.0f} | "
              f"xor_idx={r['xor_index']:.2f} sym={r['symmetry']:.2f}", flush=True)
    S = summarise(runs)
    print("\nSUMMARY (mode=%s, n=%d, shared_frac=%.2f, balanced=%d)"
          % (a.mode, S["n"], a.hub_shared_frac, a.balanced))
    print(f"  XOR fraction        : {S['xor_fraction']:.2f}")
    print(f"  both-grown fraction : {S['both_grown_fraction']:.2f}")
    print(f"  asymmetric fraction : {S['asymmetric_fraction']:.2f}")
    print(f"  R10 median [IQR]    : {S['R10'][0]:.0f} [{S['R10'][1]:.0f}, {S['R10'][2]:.0f}]")
    print(f"  R01 median [IQR]    : {S['R01'][0]:.0f} [{S['R01'][1]:.0f}, {S['R01'][2]:.0f}]")
    print(f"  R11 median [IQR]    : {S['R11'][0]:.0f} [{S['R11'][1]:.0f}, {S['R11'][2]:.0f}]")
    print(f"  xor_index median    : {S['xor_index'][0]:.2f} [{S['xor_index'][1]:.2f}, {S['xor_index'][2]:.2f}]")
    print(f"  symmetry median     : {S['symmetry'][0]:.2f} [{S['symmetry'][1]:.2f}, {S['symmetry'][2]:.2f}]")
    if a.out:
        json.dump({"summary": {k: (v if not isinstance(v, tuple) else list(v))
                               for k, v in S.items()},
                   "per_seed": [{k: (float(v) if np.isscalar(v) else None)
                                 for k, v in r.items()
                                 if k in ("seed", "R10", "R01", "R11", "R00", "A10", "A01",
                                          "xor_index", "symmetry")} for r in runs]},
                  open(a.out, "w"), indent=1)
        print("wrote", a.out)
