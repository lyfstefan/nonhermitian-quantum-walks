"""Numerical counterparts of Figs. 2-4, using finite periodic rings."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qwalk import Walk
from validation.run import match
from validation.paper_anchors import TABLES


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("examples/output"))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(3,4,figsize=(14,10),constrained_layout=True)
    for table,(model,gamma,anchors) in TABLES.items():
        row={"I":0,"II":0,"III":1,"IV":1,"V":2,"VI":2}[table]
        col=0 if model=="U1" else 2
        w=Walk(model,gamma,80)
        values,profiles,_=match(w,np.asarray(anchors))
        use=abs(w.sites)<=12
        for value,profile in zip(values,profiles.T):
            axes[row,col].plot(w.sites[use],profile[use],".-",label=f"|lambda|={abs(value):.4f}")
        axes[row,col].set(title=f"{model}, gamma={gamma}: eigenstates",xlabel="site n",ylabel="probability")
        axes[row,col].legend(fontsize=7)
        large=Walk(model,gamma,650)
        p=large.evolve()
        select=abs(large.sites)<=30
        im=axes[row,col+1].imshow(p[:,select],origin="lower",aspect="auto",extent=[-30.5,30.5,-.5,300.5],cmap="magma",vmin=0,vmax=1)
        axes[row,col+1].set(title=f"{model}, gamma={gamma}: evolution",xlabel="site n",ylabel="full step T")
        np.savez_compressed(args.output/f"{model}_gamma_{gamma}.npz",sites=large.sites,probability=p,moduli=abs(values),eigen_sites=w.sites,eigen_profiles=profiles)
    fig.colorbar(im,ax=axes[:,[1,3]],label="normalized probability",shrink=.7)
    fig.suptitle("PRA 102, 022218: 2026 numerical reimplementation (not original figure data)")
    fig.savefig(args.output/"localization.png",dpi=160)
    plt.close(fig)
    print(args.output/"localization.png")


if __name__=="__main__": main()
