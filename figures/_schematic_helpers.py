import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.lines as ml
import numpy as np
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
C_IN="#555555"; C_LAT="#1f4e79"; C_EXC="#2e7d32"; C_INH="#a11a1a"; C_DEP="#e08a00"; C_OUT="#5a3a8a"
def node(ax,xy,label,color,r=0.42,fs=11,tc="white"):
    ax.add_patch(Circle(xy,r,fc=color,ec="k",lw=1.3,zorder=5))
    ax.text(xy[0],xy[1],label,ha="center",va="center",color=tc,fontsize=fs,weight="bold",zorder=6)
    return xy,r
def ep(a,b,ra,rb):
    (x1,y1),(x2,y2)=a,b; d=np.hypot(x2-x1,y2-y1); ux,uy=(x2-x1)/d,(y2-y1)/d
    return (x1+ux*ra,y1+uy*ra),(x2-ux*rb,y2-uy*rb)
def exc(ax,A,B,color=C_EXC,lw=2.2,ls="-",ra=None,rb=None):
    p,q=ep(A[0],B[0],ra or A[1],rb or B[1])
    ax.annotate("",xy=q,xytext=p,zorder=3,arrowprops=dict(arrowstyle="-|>",color=color,lw=lw,ls=ls,shrinkA=0,shrinkB=0,mutation_scale=15))
def inh(ax,A,B,color=C_INH,lw=2.4,ra=None,rb=None,divide=False):
    p,q=ep(A[0],B[0],ra or A[1],rb or B[1])
    ax.annotate("",xy=q,xytext=p,zorder=3,arrowprops=dict(arrowstyle="-",color=color,lw=lw,shrinkA=0,shrinkB=0))
    ax.add_patch(Circle(q,0.12,fc=color,ec=color,zorder=4))
    if divide:
        mx,my=(p[0]+q[0])/2,(p[1]+q[1])/2
        ax.text(mx,my,"\u00f7",ha="center",va="center",fontsize=14,color=color,weight="bold",
                bbox=dict(boxstyle="circle,pad=0.05",fc="white",ec=color,lw=1.2),zorder=7)
