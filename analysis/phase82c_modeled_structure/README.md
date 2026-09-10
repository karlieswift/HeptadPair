# HeptadPair Phase-8.2C modeled structural case v1.0.0

Purpose: after Phase-8.2B found no exact experimental PDB structure for the selected Section-3.4 cases, generate a **qualitative modeled structural illustration** of the exact selected sequences without changing the frozen HeptadPair model or case-selection rule.

## Fixed main case

- Intended: 4H3651-AA / 4H895-AA
- Strongest competitor sharing one member: 4H3651-AA / 4H3365-AA
- Frozen XGB-1250 final shared 2434 margin: +0.156695
- Exact-recovered official set, AP=1.0

The shared member 4H3651-AA is placed as chain A in both ColabFold inputs.

## Modeling policy

Primary modeling is AlphaFold2-multimer-v3 through ColabFold with `single_sequence`, no template flag, multiple seeds and all five multimer models. This is appropriate as a qualitative structural sensitivity analysis for short designed peptides, but it is not experimental validation.

Default workload: two 29+29-residue dimers, 5 seeds x 5 models each.

## Run

```bash
cd ~/BioData/hgq/Heptad/HeptadPair_v0_1
unzip -o HeptadPair_Phase82C_ModeledStructure_v1_0_0_20260902.zip
chmod +x phase82c_modeled_structure_pack/*.sh

python phase82c_modeled_structure_pack/selftest_phase82c_analysis.py

bash phase82c_modeled_structure_pack/run_phase82c0_prepare_and_backend_preflight.sh "$PWD" \
  2>&1 | tee run_phase82c0_prepare_and_backend_preflight.log
```

If `BACKEND_PREFLIGHT.json` says `COLABFOLD_READY`:

```bash
PHASE82C_GPU=0 \
PHASE82C_NUM_SEEDS=5 \
PHASE82C_NUM_MODELS=5 \
bash phase82c_modeled_structure_pack/run_phase82c1_colabfold_main.sh "$PWD" \
  2>&1 | tee run_phase82c1_colabfold_main.log

bash phase82c_modeled_structure_pack/run_phase82c2_analyze_models.sh "$PWD" \
  2>&1 | tee run_phase82c2_analyze_models.log

bash phase82c_modeled_structure_pack/run_phase82c3_make_pymol.sh "$PWD" \
  2>&1 | tee run_phase82c3_make_pymol.log
```

If ColabFold is not in the current environment, **do not install it into `bioh200` blindly**. Use or create a separate ColabFold environment, then rerun Phase82C0. The frozen input files remain under `results/phase82c_modeled_structure/inputs/`.

## Outputs to review before manuscript figure

- `MODEL_INPUT_CONTRACT.json`
- `heptad_phase_audit.csv`
- `canonical_register_pair_audit.csv`
- `model_geometry_by_prediction.csv`
- `model_ensemble_summary.csv`
- `interchain_contacts_by_prediction.csv`
- `selected_rank1_models.csv`
- `MODELED_STRUCTURE_GATE.json`
- `pymol/*.pml` and, if PyMOL exists, rendered PNGs

## Gate

A positive geometry gate permits **qualitative modeled-structure annotation only**. It does not convert the case into an experimental structure.

## v1.0.2 hotfix (2026-09-03)

No scientific/modeling contract changes. Two post-processing bugs are fixed:

1. Score JSON lookup is now query-specific (and exact-name first). In v1.0.1,
   target and competitor files with the same `rank_###` could be cross-matched,
   contaminating reported ipTM/pTM values while leaving PDB-derived geometry,
   orientation and contact counts unchanged.
2. PyMOL script generation now uses `sel["query"]` rather than `sel.query`,
   avoiding the pandas `DataFrame.query` method-name collision.

Reuse all existing ColabFold structures. Re-run Phase82C2 and then Phase82C3 only.
