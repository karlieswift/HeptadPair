#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd


def pml_for(query, pdb, chainA_protein, chainB_protein, hmap, outpng):
    a_core=sorted(hmap[(hmap.protein==chainA_protein)&hmap.heptad.isin(list("ad"))].position.unique())
    b_core=sorted(hmap[(hmap.protein==chainB_protein)&hmap.heptad.isin(list("ad"))].position.unique())
    a_eg=sorted(hmap[(hmap.protein==chainA_protein)&hmap.heptad.isin(list("eg"))].position.unique())
    b_eg=sorted(hmap[(hmap.protein==chainB_protein)&hmap.heptad.isin(list("eg"))].position.unique())
    s=lambda xs:"+".join(map(str,xs))
    return f'''reinitialize
load {pdb}, complex
hide everything
show cartoon, complex
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.05
color cyan, chain A
color {'green' if query=='main_target' else 'orange'}, chain B
select coreA, chain A and resi {s(a_core)}
select coreB, chain B and resi {s(b_core)}
show sticks, coreA or coreB
color yellow, coreA or coreB
select egA, chain A and resi {s(a_eg)}
select egB, chain B and resi {s(b_eg)}
show sticks, egA or egB
color red, (egA or egB) and resn ASP+GLU
color blue, (egA or egB) and resn LYS+ARG
set stick_radius, 0.16
set dash_radius, 0.08
set dash_gap, 0.25
set dash_length, 0.25
bg_color white
set ray_opaque_background, off
orient complex
zoom complex, 4
ray 2400,1800
png {outpng}, dpi=300
'''


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("project_root",type=Path); args=ap.parse_args()
    root=args.project_root.resolve(); out=root/"results/phase82c_modeled_structure"
    sel=pd.read_csv(out/"selected_rank1_models.csv"); hmap=pd.read_csv(out/"sequence_heptad_map.csv")
    scripts=out/"pymol"; scripts.mkdir(exist_ok=True)
    for r in sel.itertuples(index=False):
        pml=pml_for(r.query,Path(r.pdb).resolve(),r.chainA_protein,r.chainB_protein,hmap,(scripts/f"{r.query}.png").resolve())
        (scripts/f"{r.query}.pml").write_text(pml)
    # Combined alignment: align shared chain A of competitor to target chain A and place side-by-side without claiming experimental overlay.
    t=sel[sel["query"].eq("main_target")].iloc[0]; c=sel[sel["query"].eq("main_competitor")].iloc[0]
    combo=f'''reinitialize
load {Path(t.pdb).resolve()}, target
load {Path(c.pdb).resolve()}, competitor
align competitor and chain A, target and chain A
hide everything
show cartoon, target or competitor
color cyan, target and chain A
color green, target and chain B
color cyan, competitor and chain A
color orange, competitor and chain B
set cartoon_transparency, 0.12
# Shift competitor for side-by-side comparison after alignment.
translate [35,0,0], competitor
bg_color white
set ray_opaque_background, off
orient target or competitor
zoom target or competitor, 3
ray 2800,1600
png {(scripts/'main_target_vs_competitor_side_by_side.png').resolve()}, dpi=300
'''
    (scripts/"main_target_vs_competitor_side_by_side.pml").write_text(combo)
    print(f"[OK] wrote {scripts}")

if __name__=="__main__": main()
