"""
ENGRAM-X / Brian2 proof-of-concept: the GABA switch as a moving reversal potential.
Same conductance-based E/I microcircuit, three chloride set-points:
  immature (E_Cl depolarising)  -> GDP-like synchronised population bursts (builder)
  mature   (E_Cl shunting)      -> sparse, desynchronised, stable activity (governor)
  developmental ramp            -> bursts early, vanishing at the 'GABA switch'
The only thing that changes across conditions is E_Cl - the chloride set-point.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from brian2 import *
prefs.codegen.target="numpy"; defaultclock.dt=0.5*ms

N=200; Ne=160
C=200*pF; gL=10*nS; EL=-65*mV; Ee=0*mV; Ek=-90*mV
Vth=-50*mV; Vr=-60*mV; taue=5*ms; taui=10*ms; taua=120*ms; da=4*nS
we=6*nS; wi=12*nS; BG=11*Hz
eqs='''
dv/dt=(gL*(EL-v)+ge*(Ee-v)+gi*(ECl-v)+ga*(Ek-v))/C : volt (unless refractory)
dge/dt=-ge/taue : siemens
dgi/dt=-gi/taui : siemens
dga/dt=-ga/taua : siemens
dECl/dt=rate_ecl : volt
rate_ecl : volt/second
'''
def run_cond(ecl0, ecl1, dur=4*second, seed_=2):
    start_scope(); seed(seed_)
    G=NeuronGroup(N,eqs,threshold='v>Vth',reset='v=Vr; ga+=da',refractory=3*ms,method='euler')
    G.v='EL+5*mV*rand()'; G.ECl=ecl0; G.rate_ecl=(ecl1-ecl0)/dur
    Ce=Synapses(G[:Ne],G,on_pre='ge+=we'); Ce.connect(p=0.08)
    Ci=Synapses(G[Ne:],G,on_pre='gi+=wi'); Ci.connect(p=0.10)
    P=PoissonInput(G,'ge',1,BG,weight=we)
    spk=SpikeMonitor(G); rate=PopulationRateMonitor(G); ecl=StateMonitor(G,'ECl',record=[0],dt=5*ms)
    run(dur)
    return spk, rate, ecl
def peak(rate): return float((rate.smooth_rate(window='flat',width=20*ms)/Hz).max())

sp_im,r_im,_  = run_cond(-45*mV,-45*mV)
sp_ma,r_ma,_  = run_cond(-75*mV,-75*mV)
sp_rm,r_rm,ec = run_cond(-45*mV,-75*mV)
print(f"peak population rate:  immature={peak(r_im):.0f} Hz   mature={peak(r_ma):.0f} Hz")
print(f"mean firing rate:      immature={sp_im.num_spikes/(N*4):.1f} Hz   mature={sp_ma.num_spikes/(N*4):.1f} Hz")

fig,ax=plt.subplots(2,2,figsize=(11,6.4))
ax[0,0].plot(sp_im.t/second,sp_im.i,'.',ms=1.5,color="#a11a1a")
ax[0,0].set_title("A  Immature: depolarising GABA \u2192 GDP-like bursts (builder)",weight="bold",loc="left",fontsize=9.5)
ax[0,1].plot(sp_ma.t/second,sp_ma.i,'.',ms=1.5,color="#1f4e79")
ax[0,1].set_title("B  Mature: shunting GABA \u2192 sparse & stable (governor)",weight="bold",loc="left",fontsize=9.5)
for a in [ax[0,0],ax[0,1]]: a.set_ylabel("neuron"); a.set_xlabel("time (s)"); a.set_ylim(0,N)
ax[1,0].plot(sp_rm.t/second,sp_rm.i,'.',ms=1.3,color="#555"); ax[1,0].set_ylim(0,N)
axb=ax[1,0].twinx(); axb.plot(ec.t/second,ec.ECl[0]/mV,color="#2e7d32",lw=2.2)
axb.set_ylabel("E$_{Cl}$ (mV)",color="#2e7d32"); axb.axhline(Vth/mV,ls=":",color="#999",lw=1)
axb.text(2.9,Vth/mV+1,"threshold",fontsize=7,color="#999")
ax[1,0].set_title("C  Developmental ramp: bursts vanish at the switch",weight="bold",loc="left",fontsize=9.5)
ax[1,0].set_ylabel("neuron"); ax[1,0].set_xlabel("time (s)")
mr_im=sp_im.num_spikes/(N*4); mr_ma=sp_ma.num_spikes/(N*4)
ax[1,1].bar(["immature\n(depolarising)","mature\n(shunting)"],[mr_im,mr_ma],
            color=["#a11a1a","#1f4e79"],edgecolor="k",lw=0.6)
ax[1,1].set_ylabel("mean firing rate (Hz)")
ax[1,1].set_title("D  The switch stabilises the network",weight="bold",loc="left",fontsize=9.5)
for i,v in enumerate([mr_im,mr_ma]): ax[1,1].text(i,v+0.8,f"{v:.1f} Hz",ha="center",fontsize=9)
for a in ax.ravel():
    for s in ["top","right"]: a.spines[s].set_visible(False)
fig.suptitle("ENGRAM-X in Brian2: one moving reversal potential (E$_{Cl}$) turns a builder into a governor",
             fontsize=11.5,weight="bold")
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig("demo_brian2_en.png",dpi=150); print("saved demo_brian2.png")
