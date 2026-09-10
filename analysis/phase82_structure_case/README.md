# HeptadPair Phase-8.2 structural case hotfix v1.2.0

## Why v1.1.0 stopped

The v1.1.0 generic schema detector did not recognize the real identifiers `seq1_id/seq2_id`, so the correct Phase-7.3 top-K table scored below its selection threshold. Phase-8.2B therefore never had a `selected_structure_cases.csv` to screen.

v1.2.0 removes that fragile source discovery for the primary path. It reconstructs the case directly from the frozen Phase-7.4B Section-3.4 retrieval outputs:

- `results/phase74b/locked_confirmation/family_stability/set_retrieval_metrics.csv`
- `results/phase74b/locked_confirmation/family_stability/set_partner_retrieval.csv`
- `results/phase74b/locked_confirmation/family_stability/family_repeat_averaged_predictions.parquet`

No training or tuning is performed.

## Default model

The default preference is `xgb_x2_1250_final_shared_2434` when it exists in all three frozen files, matching the strongest locked shared XGB retrieval/quantitative model. If absent, the script falls back in a declared order and records the exact chosen variant in `STRUCTURE_CASE_PREFLIGHT.json`.

## Case definition

For every official orthogonal set, intended pairs are reconstructed from the partner-retrieval table. All candidate-pair scores are taken from the already-written locked prediction table. The script independently reconstructs exact top-K recovery and audits it against the reported set metric.

- `MAIN_success_tight_competitor`: exact-recovered heterodimer with the smallest positive margin to the strongest off-target sharing either target member.
- `BACKUP_success_clear_margin`: exact-recovered heterodimer with the largest positive margin.
- `OPTIONAL_failure_boundary`: non-exact or non-positive-margin heterodimer with the smallest margin.

This is a retrospective visualization case, not a new benchmark.

## RCSB v2 hotfix

v1.2.0 also fixes the RCSB Search API v2 sequence payload before it becomes the next failure point. The query now uses:

```text
service = sequence
sequence_type = protein
identity_cutoff = 1.0
evalue_cutoff = 0.1
return_type = polymer_entity
results_content_type = experimental
```

The script then checks the RCSB Data API polymer-entity sequence. A common PDB entry is only a candidate; biological-assembly and direct-interface inspection remains mandatory.

## Run

```bash
cd ~/BioData/hgq/Heptad/HeptadPair_v0_1
unzip -o HeptadPair_Phase82_StructureCase_SelectAndPDBScreen_v1_2_0_20260902.zip
chmod +x phase82_structure_case_pack/*.sh

bash phase82_structure_case_pack/run_phase82a_structure_case_preflight.sh "$PWD" \
  2>&1 | tee run_phase82a_structure_case_preflight_v1_2_0.log

bash phase82_structure_case_pack/run_phase82b_rcsb_exact_sequence_screen.sh "$PWD" \
  2>&1 | tee run_phase82b_rcsb_exact_sequence_screen_v1_2_0.log

bash phase82_structure_case_pack/pack_phase82_results.sh "$PWD"
```

Return the resulting `HeptadPair_Phase82_StructureCase_RESULTS_v1_2_0_*.tar.gz` plus `.sha256`.
