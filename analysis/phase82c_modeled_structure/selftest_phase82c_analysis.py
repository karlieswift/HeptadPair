#!/usr/bin/env python3
from pathlib import Path
import json, math, sys, tempfile
import numpy as np
import pandas as pd
PACK=Path(__file__).resolve().parent
sys.path.insert(0,str(PACK))
from phase82c_prepare_inputs import infer_a_start
from phase82c_analyze_models import json_safe, find_score_json

seqs={
 "A":"DKNAALKAEIAALKAEIAALEAENAALEG",
 "T":"DEIAALKAKIAALKAKNAALKAKNAALEG",
 "C":"DKNAALEAKIAALEAKIAALEAKIAALEG",
}
best,_=infer_a_start(seqs)
assert best["a_start_1based"]==3, best
assert best["margin_over_second"]>5, best

probe={
    "n": np.int64(25),
    "x": np.float64(0.811),
    "arr": [np.int32(5), np.float32(0.5)],
    "nan": np.float64(np.nan),
    "path": Path('/tmp/model.pdb'),
}
safe=json_safe(probe)
assert safe["n"] == 25 and type(safe["n"]) is int, safe
assert isinstance(safe["x"], float), safe
assert safe["nan"] is None, safe
json.dumps(safe, allow_nan=False)

# Regression test: target and competitor can have the same rank number in one
# directory. The score JSON lookup must never cross query identities.
with tempfile.TemporaryDirectory() as td:
    d=Path(td)
    tp=d/'main_target__A__T_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_7.pdb'
    cp=d/'main_competitor__A__C_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_7.pdb'
    tj=d/'main_target__A__T_scores_rank_001_alphafold2_multimer_v3_model_4_seed_7.json'
    cj=d/'main_competitor__A__C_scores_rank_001_alphafold2_multimer_v3_model_1_seed_7.json'
    tp.write_text(''); cp.write_text('')
    tj.write_text('{"iptm":0.91}'); cj.write_text('{"iptm":0.72}')
    assert find_score_json(tp) == tj, (find_score_json(tp), tj)
    assert find_score_json(cp) == cj, (find_score_json(cp), cj)

print("[PASS] Phase82C heptad-phase + JSON + per-query score-JSON self-test")
