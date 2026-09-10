#!/usr/bin/env python3
import importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('p82b', HERE/'phase82b_rcsb_exact_sequence_screen.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
p=m.sequence_query_payload('A'*30)
q=p['query']['parameters']
assert q['sequence_type']=='protein'
assert q['identity_cutoff']==1.0
assert q['evalue_cutoff']==0.1
assert 'target' not in q
assert p['return_type']=='polymer_entity'
assert p['request_options']['results_content_type']==['experimental']
print('[PASS] Phase82B RCSB v2 payload self-test')
