from _helpers import *
fig=plt.figure(figsize=(9.4,8.0))
c=fig.add_axes([0.0,0.0,1.0,1.0]); c.set_xlim(0,11.6); c.set_ylim(0,10); c.axis("off")
c.set_aspect('equal')            # <-- circles stay circular (no horizontal stretch)
I2=node(c,(2.1,8.6),"I",C_IN); P2=node(c,(7.9,8.6),"P",C_IN)
Li2=node(c,(2.1,5.9),"L$_i$",C_LAT); Lp2=node(c,(7.9,5.9),"L$_p$",C_LAT)
H2=node(c,(5.0,7.1),"INH",C_INH,r=0.5)
X=node(c,(5.0,4.1),"XOR",C_OUT,r=0.54,fs=10)
Fp=node(c,(3.4,2.0),"F$^{+}$",C_EXC,r=0.4,fs=10); Fm=node(c,(6.6,2.0),"F$^{-}$",C_INH,r=0.4,fs=10)
exc(c,I2,Li2,color=C_EXC); exc(c,P2,Lp2,color=C_EXC)
exc(c,I2,H2,color=C_IN); exc(c,P2,H2,color=C_IN)
exc(c,Li2,X,color=C_EXC); exc(c,Lp2,X,color=C_EXC)
inh(c,H2,Li2,divide=True); inh(c,H2,Lp2,divide=True)
exc(c,X,Fp,color=C_EXC,lw=1.9,rb=0.4); exc(c,X,Fm,color=C_EXC,lw=1.9,rb=0.4)
exc(c,Li2,Fp,color=C_EXC,lw=1.5,rb=0.4); exc(c,Lp2,Fm,color=C_EXC,lw=1.5,rb=0.4)
c.annotate("",xy=(4.55,1.15),xytext=(3.6,1.65),arrowprops=dict(arrowstyle="-|>",color="#333",lw=1.4,mutation_scale=11))
c.annotate("",xy=(5.45,1.15),xytext=(6.4,1.65),arrowprops=dict(arrowstyle="-|>",color="#333",lw=1.4,mutation_scale=11))
c.text(5.0,0.75,"signed error  $e=\\rho(F^{+})-\\rho(F^{-})$",ha="center",fontsize=9,color="#333",
       bbox=dict(boxstyle="round,pad=0.28",fc="#f4f4f8",ec="#999"))
# divider + equations/operating cases on the right (compact)
c.plot([9.15,9.15],[1.2,9.0],color="#ddd",lw=1)
c.text(10.35,8.0,"one input on:\nlateral escapes\nthe shunt\n$\\to$ read-out ON",fontsize=7.6,ha="center",color=C_EXC)
c.text(10.35,5.3,r"$y_j=\dfrac{D_j}{\sigma+\beta\sum_k D_k}$",fontsize=11,ha="center",va="center",
       bbox=dict(boxstyle="round,pad=0.3",fc="#f4f4f8",ec="#999"))
c.text(10.35,2.7,"both on: INH $\\times2$\n$\\to$ denominator\ndivides both\n$\\to$ read-out OFF",fontsize=7.6,ha="center",color=C_INH)
fig.suptitle("The mature signed-XOR motif: shunting divides the laterals;\nthe read-out drives an opponent, signed-error pair",
             fontsize=11,weight="bold",y=0.99)
h=[ml.Line2D([],[],color=C_EXC,lw=2.4,marker=">",ms=7,label="excitatory"),
   ml.Line2D([],[],color=C_INH,lw=2.4,marker="o",ms=7,label="inhibitory / shunt (\u00f7)")]
fig.legend(handles=h,loc="lower center",ncol=2,frameon=False,fontsize=8.4,bbox_to_anchor=(0.42,0.0))
fig.savefig("motif_mature_en.png",dpi=170,bbox_inches="tight",facecolor="white")
from PIL import Image; im=Image.open("motif_mature_en.png"); print("motif_mature",im.size,"aspect",round(im.size[0]/im.size[1],2))
im.convert("RGB").save("motif_mature.jpg","JPEG",quality=92)
