#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re, os
from pathlib import Path
import pandas as pd

# Run this with the Python bundled with/able to import PyMOL, or via: pymol -cq phase82e_render_matched_views.py -- <project_root>

def resolve_path(p: str, root: Path, fallback_dir: Path):
    x=Path(str(p))
    if x.exists(): return x
    y=fallback_dir/x.name
    if y.exists(): return y
    hits=list(root.rglob(x.name))
    if len(hits)==1: return hits[0]
    raise FileNotFoundError(f'Cannot resolve {p}')

def find_rank1_candidate(root: Path, partner: str):
    ddir=root/'results/phase82d_matched_orientation'
    geom=pd.read_csv(ddir/'matched_candidate_geometry.csv')
    rows=geom[geom['query'].astype(str).str.endswith(f'__{partner}')].sort_values('rank')
    if rows.empty: raise RuntimeError(f'No Phase82D model for {partner}')
    return resolve_path(rows.iloc[0].pdb,root,ddir/'colabfold_candidates')

def rgb(hexstr):
    h=hexstr.lstrip('#'); return tuple(int(h[i:i+2],16)/255 for i in (0,2,4))

def main():
    ap=argparse.ArgumentParser(add_help=True)
    ap.add_argument('project_root', type=Path, nargs='?')
    ap.add_argument('--width', type=int, default=int(os.environ.get('PHASE82E_WIDTH','1800')))
    ap.add_argument('--height', type=int, default=int(os.environ.get('PHASE82E_HEIGHT','1450')))
    # Under `pymol -cq script.py`, PyMOL may own sys.argv. Prefer the environment variable set by the wrapper.
    args,_unknown=ap.parse_known_args(); root_arg=os.environ.get('PHASE82E_ROOT') or (str(args.project_root) if args.project_root else None)
    if not root_arg: raise SystemExit('Set PHASE82E_ROOT=/path/to/project or pass project_root when running with Python.')
    root=Path(root_arg).resolve()
    try:
        from pymol import cmd
    except Exception as e:
        raise SystemExit('PyMOL Python module is required. Run with the bioh200 PyMOL environment or: pymol -cq phase82e_render_matched_views.py -- <project_root>') from e
    cdir=root/'results/phase82c_modeled_structure'; out=root/'results/phase82e_structure_figure'; out.mkdir(parents=True,exist_ok=True)
    contract=json.loads((cdir/'MODEL_INPUT_CONTRACT.json').read_text()); a_start=int(contract['heptad_a_start_1based'])
    sel=pd.read_csv(cdir/'selected_rank1_models.csv')
    target_row=sel[sel['query'].eq('main_target')].iloc[0]; hard_row=sel[sel['query'].eq('main_competitor')].iloc[0]
    target=resolve_path(target_row.pdb,root,cdir/'colabfold_main'); hard=resolve_path(hard_row.pdb,root,cdir/'colabfold_main')
    strict=find_rank1_candidate(root,'4H140-AA'); hardpar=find_rank1_candidate(root,'4H1802-AA')
    objects={'target':target,'hard':hard,'strict':strict,'hardpar':hardpar}
    partner_color={'target':'hp_green','hard':'hp_orange','strict':'hp_magenta','hardpar':'hp_purple'}
    partner_name={'target':'4H895-AA','hard':'4H3365-AA','strict':'4H140-AA','hardpar':'4H1802-AA'}
    cmd.reinitialize()
    cmd.set_color('hp_shared',rgb('#00A9B7')); cmd.set_color('hp_green',rgb('#18A95B')); cmd.set_color('hp_orange',rgb('#E88A19'))
    cmd.set_color('hp_magenta',rgb('#B33BA5')); cmd.set_color('hp_purple',rgb('#6D54B5')); cmd.set_color('hp_core',rgb('#D8A900'))
    for name,p in objects.items(): cmd.load(str(p),name)
    # Align all shared 4H3651-AA chain-A backbones to the target. This is the only structural superposition used.
    for name in ['hard','strict','hardpar']:
        cmd.align(f'{name} and chain A and name CA', 'target and chain A and name CA', cycles=5, transform=1)
    # Canonical camera: orient on target, but zoom to union so no aligned control is clipped. Then freeze this exact view.
    cmd.hide('everything','all'); cmd.show('cartoon','target'); cmd.orient('target'); cmd.zoom('(target or hard or strict or hardpar)',6,complete=1)
    canonical_view=cmd.get_view()
    ad=[i for i in range(1,30) if 'abcdefg'[(i-a_start)%7] in 'ad']; eg=[i for i in range(1,30) if 'abcdefg'[(i-a_start)%7] in 'eg']
    adsel='+'.join(map(str,ad)); egsel='+'.join(map(str,eg))
    cmd.bg_color('white'); cmd.set('ray_opaque_background',1); cmd.set('orthoscopic',1); cmd.set('antialias',2); cmd.set('cartoon_fancy_helices',1)
    cmd.set('stick_radius',0.16); cmd.set('label_color','black'); cmd.set('label_size',18); cmd.set('label_outline_color','white')
    manifest=[]
    for name in ['target','hard','strict','hardpar']:
        cmd.hide('everything','all'); cmd.delete('hp_*label*'); cmd.delete('core*'); cmd.delete('eg*')
        cmd.show('cartoon',name); cmd.color('hp_shared',f'{name} and chain A'); cmd.color(partner_color[name],f'{name} and chain B')
        # Interface role sticks: a/d gold; acidic e/g red; basic e/g blue.
        cmd.select('coreA',f'{name} and chain A and resi {adsel}'); cmd.select('coreB',f'{name} and chain B and resi {adsel}')
        cmd.show('sticks','coreA or coreB'); cmd.color('hp_core','coreA or coreB')
        cmd.select('egA',f'{name} and chain A and resi {egsel}'); cmd.select('egB',f'{name} and chain B and resi {egsel}')
        cmd.show('sticks','egA or egB'); cmd.color('red','(egA or egB) and resn ASP+GLU'); cmd.color('blue','(egA or egB) and resn LYS+ARG')
        # N/C labels make relative orientation explicit. Shared chain always A.
        cmd.label(f'{name} and chain A and resi 1 and name CA','"A N"'); cmd.label(f'{name} and chain A and resi 29 and name CA','"A C"')
        cmd.label(f'{name} and chain B and resi 1 and name CA','"B N"'); cmd.label(f'{name} and chain B and resi 29 and name CA','"B C"')
        cmd.set_view(canonical_view)
        png=out/f'{name}_matched_view.png'; cmd.ray(args.width,args.height); cmd.png(str(png),dpi=300,ray=0)
        manifest.append({'panel':name,'partner':partner_name[name],'pdb':str(objects[name]),'png':str(png)})
    pd.DataFrame(manifest).to_csv(out/'matched_view_render_manifest.csv',index=False)
    print(f'[OK] matched-view renders -> {out}')

# PyMOL may execute a .py input with __name__ != '__main__'.
# The wrapper sets PHASE82E_PYMOL_AUTORUN=1 so the renderer runs exactly once.
if __name__ == '__main__' or os.environ.get('PHASE82E_PYMOL_AUTORUN') == '1':
    main()
