#!/usr/bin/env python3
from pathlib import Path
import tempfile,sys
import numpy as np
PACK=Path(__file__).resolve().parent
sys.path.insert(0,str(PACK))
from phase82d_analyze_orientations import orient_of
# Minimal synthetic PDB with two approximately collinear CA chains.
def atom(serial,resi,ch,x,y,z):
    return f"ATOM  {serial:5d}  CA  ALA {ch}{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 90.00           C  \n"
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'x.pdb'; lines=[]; s=1
    for i in range(1,9): lines.append(atom(s,i,'A',0,0,float(i))); s+=1
    for i in range(1,9): lines.append(atom(s,i,'B',5,0,float(i))); s+=1
    p.write_text(''.join(lines)); o,d=orient_of(p); assert o=='parallel' and d>0.9,(o,d)
print('[PASS] Phase82D orientation self-test')
