#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('p82', HERE/'phase82a_structure_case_preflight.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
variant='xgb_x2_1250_final_shared_2434'
metrics=pd.DataFrame([
 {'variant':variant,'official_set_id':'S1','exact_topk_recovery':1,'average_precision':1.0,'partner_top1_accuracy':1.0},
 {'variant':variant,'official_set_id':'S2','exact_topk_recovery':0,'average_precision':0.7,'partner_top1_accuracy':0.5},
])
partners=pd.DataFrame([
 {'variant':variant,'official_set_id':'S1','protein':'A','intended_partner':'B'},
 {'variant':variant,'official_set_id':'S1','protein':'B','intended_partner':'A'},
 {'variant':variant,'official_set_id':'S1','protein':'C','intended_partner':'D'},
 {'variant':variant,'official_set_id':'S1','protein':'D','intended_partner':'C'},
 {'variant':variant,'official_set_id':'S2','protein':'E','intended_partner':'F'},
 {'variant':variant,'official_set_id':'S2','protein':'F','intended_partner':'E'},
 {'variant':variant,'official_set_id':'S2','protein':'G','intended_partner':'H'},
 {'variant':variant,'official_set_id':'S2','protein':'H','intended_partner':'G'},
])
preds=pd.DataFrame([
 {'variant':variant,'seq_low_id':'A','seq_high_id':'B','y_pred':9},
 {'variant':variant,'seq_low_id':'C','seq_high_id':'D','y_pred':8},
 {'variant':variant,'seq_low_id':'A','seq_high_id':'C','y_pred':7.8},
 {'variant':variant,'seq_low_id':'A','seq_high_id':'D','y_pred':1},
 {'variant':variant,'seq_low_id':'B','seq_high_id':'C','y_pred':2},
 {'variant':variant,'seq_low_id':'B','seq_high_id':'D','y_pred':3},
 {'variant':variant,'seq_low_id':'E','seq_high_id':'F','y_pred':5},
 {'variant':variant,'seq_low_id':'G','seq_high_id':'H','y_pred':4},
 {'variant':variant,'seq_low_id':'E','seq_high_id':'G','y_pred':6},
 {'variant':variant,'seq_low_id':'E','seq_high_id':'H','y_pred':1},
 {'variant':variant,'seq_low_id':'F','seq_high_id':'G','y_pred':2},
 {'variant':variant,'seq_low_id':'F','seq_high_id':'H','y_pred':3},
])
_, sets, cases, _ = m.reconstruct_cases(metrics,partners,preds,variant)
assert sets.set_index('official_set_id').loc['S1','derived_exact_topk_recovery']
assert not sets.set_index('official_set_id').loc['S2','derived_exact_topk_recovery']
sel=m.choose_cases(cases)
assert 'MAIN_success_tight_competitor' in set(sel.case_role)
assert 'OPTIONAL_failure_boundary' in set(sel.case_role)
main=sel[sel.case_role.eq('MAIN_success_tight_competitor')].iloc[0]
assert abs(main.margin-0.2)<1e-8, main
print('[PASS] Phase82A selector synthetic self-test')
