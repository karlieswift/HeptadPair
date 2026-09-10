# HeptadPair Phase82E v1.0.1 - matched-view structural figure

This stage **does not retrain HeptadPair and does not rerun AlphaFold**. It reuses the frozen Phase82A-D outputs.

## Scientific roles

1. **Primary hard case**: 4H3651 + 4H895 (intended target) versus 4H3651 + 4H3365 (strongest shared-member off-target under the frozen final E+O+L XGB model). Keep this comparison primary even though the modeled orientations differ.
2. **Strict matched-orientation control**: 4H3651 + 4H140, selected only as a secondary sensitivity control because all 25 AF2-Multimer models are parallel. It does not replace the strongest off-target.
3. **Hard parallel-favored control**: 4H3651 + 4H1802, a closer frozen-score off-target whose top five models are parallel and 18/25 models are parallel. Supplementary sensitivity only.

## Why the old side-by-side was visually misleading

The previous script aligned chain A, translated the competitor, and then used a global `orient target or competitor`. The shared-chain superposition itself was good, but the single global camera, transparent background, lack of N/C labels, and occlusion made the parallel/antiparallel distinction difficult to read.

Phase82E instead aligns every candidate to the same shared chain A, freezes **one canonical camera**, renders every complex separately with exactly that view, adds N/C labels, and only then composes panels.

## Run

```bash
conda activate bioh200
cd ~/BioData/hgq/Heptad/HeptadPair_v0_1
unzip -o HeptadPair_Phase82E_StructureFigureMatchedView_v1_0_1_20260903.zip
chmod +x phase82e_structure_figure_pack/*.sh phase82e_structure_figure_pack/*.py
bash phase82e_structure_figure_pack/run_phase82e_structure_figure.sh "$PWD" \
  2>&1 | tee run_phase82e_structure_figure.log
```

Outputs are written to `results/phase82e_structure_figure/`.

## Wording boundary

These are **modeled/predicted structures**, not experimental structures. Model/seed ensemble members are descriptive structural samples, not independent biological replicates. Do not use them for inferential p-values.


## v1.0.1 hotfix

The v1.0.0 wrapper called `pymol -cq phase82e_render_matched_views.py`, but PyMOL can execute an input Python file without setting `__name__ == "__main__"`. The renderer therefore could be loaded without calling `main()`, after which the composer failed because `target_matched_view.png` did not exist.

v1.0.1 sets `PHASE82E_PYMOL_AUTORUN=1`, explicitly autoruns the renderer under PyMOL, and adds a hard output gate requiring all four matched-view PNGs plus the render manifest before composition. No scientific selection, PDB coordinates, metrics, or figure interpretation were changed. No AlphaFold rerun is required.
