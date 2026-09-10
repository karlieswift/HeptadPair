#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from PIL import Image, ImageChops

ap=argparse.ArgumentParser()
ap.add_argument('--repo',default='.')
ap.add_argument('--outdir',default='reproduced/figures')
args=ap.parse_args()
ROOT=Path(args.repo).resolve()
out_arg=Path(args.outdir)
OUT=out_arg if out_arg.is_absolute() else ROOT/out_arg
MAIN=OUT; SUP=OUT
MAIN.mkdir(parents=True,exist_ok=True); SUP.mkdir(parents=True,exist_ok=True)

# Palette consistent with existing manuscript figures
blue='#2E73B8'; teal='#1F9E89'; purple='#7A4EB2'; orange='#E98920'; navy='#143F69'; gray='#697386'; lightgray='#D8DEE8'; red='#D84A3A'; green='#239E84'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10})

def save(fig,pdf,png,dpi=300):
    fig.savefig(pdf,bbox_inches='tight')
    fig.savefig(png,dpi=dpi,bbox_inches='tight')
    plt.close(fig)

def add_panel(ax,letter,x=-0.08,y=1.06):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=15,fontweight='bold',va='top',ha='left',color='#20242A')

def box(ax,xy,w,h,text,edge,color_text=None,sub=None):
    p=FancyBboxPatch(xy,w,h,boxstyle='round,pad=0.012,rounding_size=0.015',fc='white',ec=edge,lw=1.7)
    ax.add_patch(p)
    ax.text(xy[0]+w/2,xy[1]+h*0.58,text,ha='center',va='center',fontweight='bold',fontsize=9.5,color=color_text or edge)
    if sub:
        ax.text(xy[0]+w/2,xy[1]-0.04,sub,ha='center',va='top',fontsize=7.6,color=gray)
    return p

# ---------------- Figure 2 ----------------
fact=pd.read_csv(ROOT/'results/phase81_factorial_ablation_v1_0_2/factorial_metrics_by_split.csv')
fsum=pd.read_csv(ROOT/'results/phase81_factorial_ablation_v1_0_2/factorial_summary.csv')
loo=pd.read_csv(ROOT/'results/phase81_factorial_ablation_v1_0_2/leave_one_block_out_effects_by_split.csv')
loos=pd.read_csv(ROOT/'results/phase81_factorial_ablation_v1_0_2/leave_one_block_out_summary.csv')

fig=plt.figure(figsize=(14,10.6))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.03],width_ratios=[1.08,0.92],hspace=0.30,wspace=0.16)
fig.suptitle('Figure 2 | From seven-center graphs to the frozen 2,434-dimensional HeptadPair representation',x=0.05,y=0.985,ha='left',fontsize=18,fontweight='bold')

# a development path
ax=fig.add_subplot(gs[0,0]); ax.set_axis_off(); add_panel(ax,'a',-0.03,1.03)
ax.text(0.04,0.98,'Biological compression, not graph complexity,\ndrove model development',transform=ax.transAxes,ha='left',va='top',fontsize=11.0,fontweight='bold',color=navy,linespacing=1.05)
ax.set_xlim(0,1); ax.set_ylim(0,1)
coords={
 'g':(0.06,0.55), 't':(0.35,0.55), 'e':(0.64,0.55),
 'o':(0.64,0.18), 'l':(0.35,0.18), 'x':(0.06,0.18)}
box(ax,coords['g'],.22,.17,'1  Seven-center\nGNN',gray,gray,'all seven hubs')
box(ax,coords['t'],.22,.17,'2  Typed-four\ntwo-center',blue,blue,'a/d core + e/g charge')
box(ax,coords['e'],.22,.17,'3  Engineered\n2,214',teal,teal,'deterministic base')
box(ax,coords['o'],.22,.17,'4  Ordered-register\n+118',purple,purple,'global typed order')
box(ax,coords['l'],.22,.17,'5  Local-register\n+102',purple,purple,'local register context')
box(ax,coords['x'],.22,.17,'6  XGB-1250\nlocked CCmax',navy,navy,'task-specific calibration')

def arrow(a,b):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,lw=1.2,color='#2C3640',shrinkA=2,shrinkB=2))
arrow((.28,.635),(.35,.635)); arrow((.57,.635),(.64,.635)); arrow((.75,.55),(.75,.35)); arrow((.64,.265),(.57,.265)); arrow((.35,.265),(.28,.265))
ax.text(.04,.04,'7-center graph → typed-four graph → 2,214 + 118 + 102 = 2,434 → locked learner',fontsize=9.2,color='#223044')

# b frozen benchmark
ax=fig.add_subplot(gs[0,1]); add_panel(ax,'b',-0.07,1.05)
ax.set_title('Frozen CCNG1 both-chains-unseen benchmark',loc='left',fontweight='bold',color=navy,pad=7)
labels=['iCipa\nreimpl.','Engineered\nHGB 2,214','Ordered\nHGB 2,332','Final\nHGB 2,434']
pear=[0.2798,0.6352,0.6572,0.6639]; spear=[0.2502,0.6014,0.6263,0.6305]
x=np.arange(len(labels)); w=.34
ax.bar(x-w/2,pear,w,label='Pearson',color=blue,alpha=.96)
ax.bar(x+w/2,spear,w,label='Spearman',color=teal,alpha=.96)
for i,v in enumerate(spear): ax.text(i+w/2,v+0.008,f'{v:.3f}',ha='center',fontsize=8.6)
ax.set_xticks(x,labels); ax.set_ylim(0,0.72); ax.set_ylabel('Correlation'); ax.grid(axis='y',alpha=.18); ax.legend(frameon=False,loc='upper left')
for s in ['top','right']: ax.spines[s].set_visible(False)

# c factorial
ax=fig.add_subplot(gs[1,0]); add_panel(ax,'c',-0.05,1.07)
ax.set_title('Full-factorial feature-block ablation',loc='left',fontsize=15,fontweight='bold',color=navy,pad=8)
order=['E','O','L','EO','EL','OL','EOL']; dims={'E':2214,'O':118,'L':102,'EO':2332,'EL':2316,'OL':220,'EOL':2434}
cols={'E':blue,'O':purple,'L':orange,'EO':'#4C77A8','EL':'#4F8D98','OL':'#87618E','EOL':teal}
for i,sub in enumerate(order):
    g=fact[(fact['subset']==sub)].sort_values('split_seed')
    vals=g['spearman'].to_numpy()
    # deterministic symmetric jitter
    jit=np.linspace(-0.10,0.10,len(vals))
    ax.scatter(np.full(len(vals),i)+jit,vals,s=42,facecolors='white',edgecolors=cols[sub],linewidths=1.6,zorder=3)
    m=np.nanmean(vals); sd=np.nanstd(vals,ddof=1)
    ax.errorbar(i,m,yerr=sd,fmt='D',ms=8,color=cols[sub],capsize=4,lw=2,zorder=4)
    ax.text(i,m+0.014,f'{m:.3f}',ha='center',fontsize=8.8,color='#26313E')
ax.axvline(2.5,color=lightgray,lw=1); ax.axvline(5.5,color=lightgray,lw=1)
ax.text(1,0.472,'single blocks',ha='center',color=gray,fontsize=9)
ax.text(4,0.472,'two-block combinations',ha='center',color=gray,fontsize=9)
ax.text(6,0.472,'full',ha='center',color=gray,fontsize=9)
ax.set_xticks(range(7),[f'{s}\n{dims[s]:,}' for s in order]); ax.set_ylabel('Test Spearman correlation'); ax.set_ylim(0.46,0.725); ax.grid(axis='y',alpha=.18)
ax.text(4.0,0.495,'O+L (220) ≈ E (2,214)',ha='center',va='center',fontsize=8.8,color=navy,bbox=dict(boxstyle='round,pad=.18',fc='white',ec='#C8D6E6'))
for s in ['top','right']: ax.spines[s].set_visible(False)

# d LOO
ax=fig.add_subplot(gs[1,1]); add_panel(ax,'d',-0.08,1.07)
ax.set_title('Leave-one-block-out contribution',loc='left',fontsize=15,fontweight='bold',color=navy,pad=8)
ypos={'E':2,'O':1,'L':0}; blocknames={'E':'E: engineered','O':'O: ordered-register','L':'L: local-register'}; ccols={'E':blue,'O':purple,'L':orange}
for b in ['E','O','L']:
    g=loo[(loo.block==b)&(loo.metric=='spearman')].sort_values('split_seed')
    y=ypos[b]
    jit=np.linspace(-.12,.12,len(g))
    ax.scatter(g['oriented_gain'],np.full(len(g),y)+jit,s=40,facecolors='white',edgecolors=ccols[b],linewidths=1.6,zorder=3)
    s=loos[(loos.block==b)&(loos.metric=='spearman')].iloc[0]
    mean=s.oriented_gain_mean; lo=s.bootstrap_ci95_low; hi=s.bootstrap_ci95_high
    ax.errorbar(mean,y,xerr=[[mean-lo],[hi-mean]],fmt='D',ms=8,color=ccols[b],capsize=4,lw=2.2,zorder=4)
    ax.text(0.079,y,f'+{mean:.3f}',ha='right',va='center',fontsize=10,fontweight='bold',color=ccols[b])
ax.axvline(0,color='#20242A',ls='--',lw=1.2)
ax.set_yticks([2,1,0],['E','O','L']); ax.set_xlim(-0.045,0.085); ax.set_xlabel('ΔSpearman: full E+O+L minus reduced model'); ax.grid(axis='x',alpha=.18)
ax.text(-0.043,-0.50,'Five frozen outer splits; bars are bootstrap 95% CIs.\nSmall n makes exact sign-flip p-values coarse; effects are descriptive.',fontsize=8.4,color=gray,va='top')
for s in ['top','right','left']: ax.spines[s].set_visible(False)
ax.tick_params(axis='y',length=0)

save(fig,MAIN/'Figure2_Model_Evolution_CCNG1.pdf',MAIN/'Figure2_Model_Evolution_CCNG1.png')

# ---------------- Figure 4 ----------------
case=pd.read_csv(ROOT/'results/phase82e_structure_figure/phase82e_case_summary.csv').set_index('label')

def crop_nonwhite(path, pad=15):
    im=Image.open(path).convert('RGB')
    bg=Image.new('RGB',im.size,'white')
    diff=ImageChops.difference(im,bg).convert('L')
    bbox=diff.point(lambda p: 255 if p>18 else 0).getbbox()
    if not bbox: return im
    l,t,r,b=bbox; l=max(0,l-pad); t=max(0,t-pad); r=min(im.width,r+pad); b=min(im.height,b+pad)
    return im.crop((l,t,r,b))

structs={
 'target':crop_nonwhite(ROOT/'results/phase82e_structure_figure/target_matched_view.png'),
 'hard':crop_nonwhite(ROOT/'results/phase82e_structure_figure/hard_matched_view.png'),
 'strict':crop_nonwhite(ROOT/'results/phase82e_structure_figure/strict_matched_view.png'),
 'hardpar':crop_nonwhite(ROOT/'results/phase82e_structure_figure/hardpar_matched_view.png')
}

fig=plt.figure(figsize=(14,12.6))
gs=fig.add_gridspec(3,4,height_ratios=[1.0,0.86,1.12],hspace=0.30,wspace=0.30)
fig.suptitle('Figure 4 | Orthogonal partner recovery and a modeled hard-pair structural case',x=.055,y=.988,ha='left',fontsize=18,fontweight='bold')
# a
ax=fig.add_subplot(gs[0,0:2]); add_panel(ax,'a',-.06,1.08); ax.set_title('Engineered features lead zero-shot retrieval',loc='left',fontweight='bold',color=navy,pad=6)
metrics=['AP','Recall@K','Top-1','MRR']; eng=[.8595,.7935,.8286,.8765]; fin=[.7097,.6228,.6324,.7605]
x=np.arange(4); w=.34
ax.bar(x-w/2,eng,w,label='Engineered 2,214',color=blue); ax.bar(x+w/2,fin,w,label='Final 2,434',color=teal)
ax.set_xticks(x,metrics); ax.set_ylim(.55,.92); ax.set_ylabel('Zero-shot retrieval metric'); ax.grid(axis='y',alpha=.18); ax.legend(frameon=False,fontsize=8)
for ss in ['top','right']: ax.spines[ss].set_visible(False)
# b
ax=fig.add_subplot(gs[0,2:4]); add_panel(ax,'b',-.06,1.08); ax.set_title('Calibration improves orthogonal-set discrimination',loc='left',fontweight='bold',color=navy,pad=6)
vals=[.8595,.7097,.9568]; labs=['Zero-shot\nengineered','Zero-shot\nfinal','Supervised\ncalibrated']; cc=[blue,teal,purple]
bars=ax.bar(np.arange(3),vals,color=cc,width=.58)
for bb,v in zip(bars,vals): ax.text(bb.get_x()+bb.get_width()/2,v+.008,f'{v:.4f}',ha='center',fontweight='bold',fontsize=8.5)
ax.set_ylim(.65,1.0); ax.set_ylabel('Average precision'); ax.set_xticks(np.arange(3),labs); ax.grid(axis='y',alpha=.18)
for ss in ['top','right']: ax.spines[ss].set_visible(False)
# c
ax=fig.add_subplot(gs[1,0:2]); add_panel(ax,'c',-.06,1.08); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.text(.04,.95,'Best exact recovery: 9 of 13 predefined sets',fontsize=11,fontweight='bold',color=navy,va='top')
for i in range(13):
    row=0 if i<7 else 1; col=i if i<7 else i-7
    x0=.06+col*.125; y0=.51-row*.28
    fc=teal if i<9 else 'white'; ec=teal if i<9 else '#8590A0'; txt='white' if i<9 else '#596273'
    rr=FancyBboxPatch((x0,y0),.09,.17,boxstyle='round,pad=.008,rounding_size=.012',fc=fc,ec=ec,lw=1.3)
    ax.add_patch(rr); ax.text(x0+.045,y0+.085,str(i+1),ha='center',va='center',fontweight='bold',color=txt,fontsize=9)
ax.text(.06,.08,'Filled tiles indicate complete intended-pair recovery;\nset identities are provided in Source Data.',fontsize=8.2,color=gray)
# d
ax=fig.add_subplot(gs[1,2:4]); add_panel(ax,'d',-.06,1.08); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.text(.03,.96,'Orthogonality is a top-end retrieval endpoint',fontsize=11,fontweight='bold',color=navy,va='top')
for y,title,body,edge,fc in [
 (.51,'Biological interpretation','Intended pairs must outrank all admissible off-targets.\nThe a/d core and e/g charge channels contribute\ncomplementary constraints to this ranking.',orange,'#FFF8EF'),
 (.11,'Primary limitation','Performance decreased as homodimer fraction increased\n(r = -0.811). Self-compatibility and partner asymmetry\nremain incompletely resolved.',red,'#FFF2F1')]:
    rr=FancyBboxPatch((.05,y),.90,.28,boxstyle='round,pad=.015,rounding_size=.02',fc=fc,ec=edge,lw=1.5); ax.add_patch(rr)
    ax.text(.10,y+.22,title,fontweight='bold',color=edge,fontsize=10,va='top'); ax.text(.10,y+.14,body,fontsize=8.2,color='#343A40',va='top')
# e full-width matched structural comparison
subgs=gs[2,:].subgridspec(1,3,wspace=.06)
roles=[('target','Intended target\n4H3651 + 4H895','target'),('hard','Strongest shared-member off-target\n4H3651 + 4H3365','hard_competitor'),('strict','Strict matched-parallel control\n4H3651 + 4H140','strict_parallel_control')]
for j,(key,title,rowname) in enumerate(roles):
    ax=fig.add_subplot(subgs[0,j]);
    if j==0: add_panel(ax,'e',-.08,1.08)
    ax.imshow(structs[key]); ax.axis('off'); ax.set_title(title,fontsize=10.0,fontweight='bold',pad=2)
    r=case.loc[rowname]
    orient='parallel' if r.parallel_n==25 else ('antiparallel' if r.antiparallel_n==25 else 'mixed')
    stats=(f"E+O+L XGB {r.frozen_eol_xgb_score:.3f} | {orient}\n"
           f"mean ipTM {r.mean_iptm:.3f}; median contacts {r.median_contacts:.0f}; "
           f"a/d {r.median_ad:.0f}; e/g opp/same {r.median_eg_opp:.0f}/{r.median_eg_same:.0f}")
    ax.text(.5,-.025,stats,transform=ax.transAxes,ha='center',va='top',fontsize=7.9,color='#2F3742',linespacing=1.25)
fig.text(.055,.014,'Shared chain 4H3651 is aligned to an identical camera. AlphaFold2-Multimer/ColabFold coordinates are retrospective modeled illustrations, not experimental structures; model/seed samples are descriptive.',fontsize=8.2,color=gray,ha='left')
save(fig,MAIN/'Figure4_Orthogonal_Partner_Recovery.pdf',MAIN/'Figure4_Orthogonal_Partner_Recovery.png')

# ---------------- Supplementary Figure S6 ----------------
fig=plt.figure(figsize=(11.8,5.6))
gs=fig.add_gridspec(1,3,width_ratios=[1.25,1.25,0.9],wspace=.16)
fig.suptitle('Supplementary Figure S6 | Hard parallel-favored off-target sensitivity control',x=.05,y=.98,ha='left',fontsize=15,fontweight='bold')
for j,(key,title) in enumerate([('target','Intended target\n4H3651 + 4H895'),('hardpar','Hard parallel-favored off-target\n4H3651 + 4H1802')]):
    ax=fig.add_subplot(gs[0,j]); ax.imshow(structs[key]); ax.axis('off'); ax.set_title(title,fontsize=10,fontweight='bold')
    if j==0: ax.text(-.06,1.04,'a',transform=ax.transAxes,fontweight='bold',fontsize=14)
ax=fig.add_subplot(gs[0,2]); ax.set_axis_off(); ax.text(-.10,1.04,'b',transform=ax.transAxes,fontweight='bold',fontsize=14)
ax.text(.02,.98,'Frozen-score and ensemble summary',fontweight='bold',color=navy,fontsize=10,va='top')
# hardpar parallel subset summary
hp=pd.read_csv(ROOT/'results/phase82e_structure_figure/phase82e_4H1802_parallel_subset_summary.csv').iloc[0]
lines=[
 ('Target score','-0.511'),('4H1802 score','-0.772'),('4H1802 orientation','18/25 parallel'),('Top-5 orientation','5/5 parallel'),('Parallel-subset mean ipTM',f"{hp.mean_iptm:.3f}"),('median contacts',f"{hp.median_contacts:.0f}"),('median a/d',f"{hp.median_ad:.0f}"),('median e/g opp./same',f"{hp.median_eg_opp:.0f}/{hp.median_eg_same:.0f}")]
for i,(k,v) in enumerate(lines):
    y=.86-i*.095; ax.text(.03,y,k,fontsize=8.2,color='#404852'); ax.text(.97,y,v,fontsize=8.2,ha='right',fontweight='bold')
ax.add_patch(FancyBboxPatch((.03,.035),.94,.13,boxstyle='round,pad=.01',fc='#F5F7FA',ec='#BBC5D1'))
ax.text(.06,.14,'Boundary',fontweight='bold',fontsize=8,color=navy,va='top')
ax.text(.06,.10,'Orientation is not 25/25 stable; this case\nis therefore supplementary and does not\nreplace the strict matched-parallel control.',fontsize=7.5,va='top',color='#38404A')
fig.text(.05,.02,'Shared chain 4H3651 is aligned to the same camera. AlphaFold2-Multimer models are structural sensitivity illustrations only.',fontsize=8.2,color=gray)
save(fig,SUP/'Supplementary_Figure_S6_Phase82_Structural_Control.pdf',SUP/'Supplementary_Figure_S6_Phase82_Structural_Control.png')

print('WROTE',MAIN/'Figure2_Model_Evolution_CCNG1.pdf')
print('WROTE',MAIN/'Figure4_Orthogonal_Partner_Recovery.pdf')
print('WROTE',SUP/'Supplementary_Figure_S6_Phase82_Structural_Control.pdf')
