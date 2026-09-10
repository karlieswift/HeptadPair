# Reproducibility guide

## Tier 1 — verify the frozen release

Requires only Python 3.11+:

```bash
python scripts/repo_doctor.py --repo . --strict
python scripts/verify_release.py --repo .
```

The first command checks required files, path hygiene and release gates. The second verifies every checksum recorded in `MANIFEST.tsv`.

## Tier 2 — validate source data and rebuild the revised figures

Create the lightweight environment:

```bash
conda env create -f environment/environment.yml
conda activate heptadpair_release
python scripts/check_source_data.py --repo .
python scripts/reproduce_revision_figures.py --repo . --outdir reproduced/figures
```

This rebuilds revised Figure 2, Figure 4 and Supplementary Figure S6 from the included public-safe tables and matched-view structural renders. Figures 1, 3 and 5 are retained as frozen manuscript artifacts.

## Tier 3 — re-run Phase82 coordinate summaries from the included modeled coordinates

The compact release includes the 5-seed x 5-model PDB/score outputs for the target, strongest off-target and matched-orientation candidates. The exact Phase82 analysis scripts are under `analysis/` and the frozen outputs are under `results/`.

The original Phase82 wrapper scripts were written for the historical project layout. For public verification, use the included precomputed model-level tables and Tier-2 reconstruction unless you intentionally reproduce the historical directory contract.

## Tier 4 — regenerate AlphaFold2-Multimer coordinates

Model weights and third-party databases are not redistributed. The exact software-version snapshot is documented in `environment/COLABFOLD82C_ENVIRONMENT_SANITIZED.txt`, and the input/analysis contract is in `contracts/STRUCTURAL_CASE_CONTRACT.md`.

## Phase81 from-scratch note

The Phase81 scripts are included, but a from-scratch replay requires the historical Phase67 feature tables, exact split files and local HeptadPair package. Those large/historical project assets are not duplicated in this compact release. The release therefore provides the complete frozen Phase81 output tables plus the Phase67 replay/split-reuse audits.
