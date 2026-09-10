#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root',type=Path); args=ap.parse_args(); root=args.project_root.resolve()
    out=root/'results/phase82e_structure_figure'; s=pd.read_csv(out/'phase82e_case_summary.csv').set_index('label')
    panels=[('target','Intended target\n4H3651 + 4H895'),('hard','Strongest shared-member off-target\n4H3651 + 4H3365'),('strict','Strict matched-parallel control\n4H3651 + 4H140')]
    labelmap={'target':'target','hard':'hard_competitor','strict':'strict_parallel_control'}
    fig,axs=plt.subplots(1,3,figsize=(14.7,5.8),constrained_layout=True)
    for ax,(p,title) in zip(axs,panels):
        img=mpimg.imread(out/f'{p}_matched_view.png'); ax.imshow(img); ax.axis('off'); ax.set_title(title,fontsize=11,fontweight='bold')
        r=s.loc[labelmap[p]]
        ori=('parallel' if r.parallel_n==r.n_models else ('antiparallel' if r.antiparallel_n==r.n_models else f"{int(r.parallel_n)}/{int(r.n_models)} parallel"))
        txt=(f"E+O+L XGB score {r.frozen_eol_xgb_score:.3f}\n"
             f"AF2 orientation {ori}; mean ipTM {r.mean_iptm:.3f}\n"
             f"median contacts {r.median_contacts:.0f}; a/d {r.median_ad:.0f}; e/g opp/same {r.median_eg_opp:.0f}/{r.median_eg_same:.0f}")
        ax.text(0.5,-0.03,txt,transform=ax.transAxes,ha='center',va='top',fontsize=8.5,linespacing=1.25)
    fig.suptitle('Modeled partner-specific structural case: same query chain, distinct candidate partners',fontsize=13,fontweight='bold')
    fig.text(0.5,0.005,'All panels use the same shared-chain alignment and identical camera. AlphaFold2-Multimer models are retrospective structural illustrations, not experimental structures.',ha='center',fontsize=8.5)
    fig.savefig(out/'Figure4_structure_case_triptych.png',dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(out/'Figure4_structure_case_triptych.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)
    # Supplementary harder, orientation-bimodal control.
    fig,axs=plt.subplots(1,2,figsize=(10.2,5.7),constrained_layout=True)
    for ax,p,title,label in [(axs[0],'target','Intended target: 4H3651 + 4H895','target'),(axs[1],'hardpar','Hard parallel-favored control: 4H3651 + 4H1802','hard_parallel_favored_control')]:
        ax.imshow(mpimg.imread(out/f'{p}_matched_view.png')); ax.axis('off'); ax.set_title(title,fontsize=11,fontweight='bold')
        r=s.loc[label]
        ax.text(0.5,-0.03,f"E+O+L XGB score {r.frozen_eol_xgb_score:.3f}; orientation {int(r.parallel_n)}/{int(r.n_models)} parallel; mean ipTM {r.mean_iptm:.3f}",transform=ax.transAxes,ha='center',va='top',fontsize=8.5)
    fig.savefig(out/'Supplementary_structure_parallel_favored_control.png',dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(out/'Supplementary_structure_parallel_favored_control.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)
    print(f'[OK] composed figures -> {out}')
if __name__=='__main__': main()
