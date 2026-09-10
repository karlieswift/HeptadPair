# Revision analysis scripts

This directory preserves the exact Phase81 and Phase82 revision-stage scripts used to generate the frozen results.

- `phase81_full_factorial/` — complete E/O/L factorial ablation code. A from-scratch run requires historical Phase67 features, split files and the local HeptadPair package.
- `phase82_structure_case/` — score-first structural-case selection and exact experimental-PDB sequence screen.
- `phase82c_modeled_structure/` — modeled-structure preparation, ColabFold run wrapper, coordinate/contact analysis and PyMOL scripts.
- `phase82d_matched_orientation/` — score-first matched-orientation control selection/modeling audit.
- `phase82e_structure_figure/` — matched-view structural audit/rendering/composition.

For a compact reviewer workflow, use the portable scripts under `scripts/` against the included frozen outputs rather than reconstructing the full historical project tree.
