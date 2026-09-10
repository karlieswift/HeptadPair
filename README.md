# HeptadPair

**A cyclic register-aware representation for unseen coiled-coil interaction prediction and orthogonal partner recovery**

This repository is a compact, reviewer-facing reproducibility snapshot for the v8.9 HeptadPair manuscript revision. It integrates the frozen Phase81 full-factorial E/O/L analysis and the Phase82 modeled partner-specific structural case without reopening model selection.

![HeptadPair representation](figures/main/Figure1_Cyclic_RegisterAware_TypedFour_DualCenter.png)

## What changed in v8.9

### Full-factorial E/O/L analysis

All seven non-empty combinations of engineered (E), ordered-register (O) and local-register (L) blocks were evaluated on the same five frozen CCNG1 both-chains-unseen splits. The complete E+O+L representation gives the highest mean Spearman (0.6321). Leave-one-block-out losses support the conditional ranking **E > O > L**. The compact O+L block (220 dimensions) nearly matches the 2,214-dimensional E block.

### Partner-specific modeled structural case

The final frozen E+O+L XGB model ranks the intended pair `4H3651-AA + 4H895-AA` above the strongest shared-member off-target `4H3651-AA + 4H3365-AA`. Both form confident modeled coiled coils; the off-target actually has more modeled total and a/d contacts. A score-first strict parallel control, `4H3651-AA + 4H140-AA`, has higher modeled ipTM and denser a/d packing than the target but a substantially lower HeptadPair score. The structural models therefore illustrate **partner-specific discrimination beyond generic modeled foldability/interface size/core packing**.

No exact experimental PDB structure was found for the selected pair. All Phase82 coordinates are modeled/predicted illustrations; ipTM/pTM are confidence measures, not affinity.

## Repository map

- `manuscript/` — final v8.9 LaTeX/PDF manuscript.
- `supplementary/` — final v8.9 Supplementary Information.
- `figures/` — main Figures 1–5 and Supplementary Figures S1–S6.
- `source_data/` — frozen base-manuscript source-data workbook plus the public-safe v8.9 Phase81/Phase82 workbook and CSV exports.
- `contracts/` — frozen dataset, feature and structural-case contracts.
- `results/` — frozen Phase81/Phase82 tables, gates and modeled-coordinate outputs.
- `analysis/` — exact Phase81/Phase82 revision scripts.
- `scripts/` — portable release verification and revised-figure reconstruction.
- `environment/` — lightweight reviewer environment and exact ColabFold environment provenance.
- `docs/` — data boundary, result index, claim scope and reproducibility guide.

## Quick verification

```bash
python scripts/repo_doctor.py --repo . --strict
python scripts/verify_release.py --repo .
```

For the revised figures/source data:

```bash
conda env create -f environment/environment.yml
conda activate heptadpair_release
python scripts/check_source_data.py --repo .
python scripts/reproduce_revision_figures.py --repo . --outdir reproduced/figures
```

See `docs/REPRODUCIBILITY.md` for the full hierarchy of frozen-output verification versus historical from-scratch replay.

## Scientific boundaries

- The Z7 construction encodes heptad register/relative phase rather than exact 3D geometry.
- The current framework targets dimeric coiled-coil pair compatibility, not higher-order oligomerization.
- “Unseen” means absent from the corresponding model-training set, not newly discovered/generated de novo.
- Learner and feature subset are endpoint-specific; the framework does not claim one universal learner.
- The Phase82 structural case is a retrospective illustration of a frozen retrieval decision, not a new independent validation benchmark.
- The production 2,434-dimensional representation excludes the non-promoted GLMY/path-homology audit.

## Data and code boundary

Processed source tables, frozen revision results, revision scripts and modeled Phase82 outputs needed for audit are included. Large raw third-party datasets and the complete historical training repository are not duplicated; obtain upstream public data from the sources cited in the manuscript.

## Release status

`PUBLIC_RELEASE_GATE.json` reports:

- `reviewer_ready = true`
- `github_package_ready = true`
- `public_release_ready = false`

The remaining blockers are author-controlled release metadata, especially license selection, repository/archival identifiers, and final administrative statements. No license is inferred by this package.
