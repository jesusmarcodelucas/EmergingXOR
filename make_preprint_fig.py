import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from matplotlib.patches import FancyArrowPatch
plt.rcParams.update({"font.family":"DejaVu Serif","font.size":9})
fig,ax=plt.subplots(figsize=(9.2,3.9))
x=np.linspace(0,10,400)
df=40*np.exp(-0.55*x)                      # driving force for GABA (mV above rest): depol -> ~0
ax.plot(x,df,color="#a11a1a",lw=2.4)
ax.axhline(0,color="#333",lw=1.0)
ax.fill_between(x,0,df,where=(df>4),color="#f6dcdc",alpha=0.7)
sw=np.argmin(np.abs(df-4)); xsw=x[sw]
ax.axvline(xsw,color="#2e7d32",ls="--",lw=1.4)
# GDP frequency (secondary)
ax2=ax.twinx()
gdp=np.exp(-((x-3.2)**2)/2.2)*np.clip((xsw-x)/xsw,0,1)*1.4
ax2.plot(x,gdp,color="#1f4e79",lw=1.8,ls="-")
ax2.set_ylim(0,2.2); ax2.set_ylabel("GDP frequency (a.u.)",color="#1f4e79",fontsize=8.5)
ax2.tick_params(axis='y',labelcolor="#1f4e79",labelsize=7)
for s in ["top"]: ax2.spines[s].set_visible(False)
# labels
ax.set_xlim(0,10); ax.set_ylim(-6,46)
ax.set_xlabel("developmental age  (embryonic  \u2192  birth  \u2192  postnatal)")
ax.set_ylabel("GABA driving force (mV above rest)",color="#a11a1a",fontsize=8.5)
ax.tick_params(axis='y',labelcolor="#a11a1a"); ax.set_xticks([])
ax.text(2.6,42,"depolarising GABA  (excitatory)",fontsize=8.2,color="#a11a1a")
ax.text(7.2,3.5,"shunting GABA  (inhibitory)",fontsize=8.2,color="#7a2a2a")
ax.text(xsw+0.1,30,"GABA switch\nKCC2 \u2191 , [Cl\u207b] \u2193",fontsize=8,color="#2e7d32",va="top")
# stage annotations
ax.annotate("interneurons &\nGABA synapses first",xy=(0.6,8),xytext=(0.3,-5.3),fontsize=7.4,color="#333",ha="left")
ax.annotate("depolarising GABA \u2192 GDPs\nbuild the glutamatergic arms",xy=(3.2,1.5),xytext=(2.2,18),fontsize=7.6,color="#1f4e79",ha="left",
            arrowprops=dict(arrowstyle="->",color="#1f4e79",lw=0.9))
ax.annotate("mature shunting\nsigned-XOR",xy=(8.6,1.0),xytext=(8.0,20),fontsize=7.6,color="#7a2a2a",ha="left")
# builder / governor mode bars
ax.annotate("",xy=(xsw-0.1,44),xytext=(0.2,44),arrowprops=dict(arrowstyle="-",color="#c05a5a",lw=6,alpha=0.35))
ax.annotate("",xy=(9.9,44),xytext=(xsw+0.1,44),arrowprops=dict(arrowstyle="-",color="#4a7a4a",lw=6,alpha=0.35))
ax.text((xsw)/2,45.6,"INH hub = BUILDER (cooperative)",fontsize=7.6,ha="center",color="#8a2a2a",weight="bold")
ax.text((xsw+10)/2,45.6,"INH hub = GOVERNOR (shunting / XOR fold)",fontsize=7.6,ha="center",color="#2e6b2e",weight="bold")
for s in ["top"]: ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig("preprint_fig_en.png",dpi=180,bbox_inches="tight")
from PIL import Image; print("saved",Image.open("preprint_fig_en.png").size)
