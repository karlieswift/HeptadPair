# HeptadPair Phase81 full-factorial ablation v1.0.2

Hotfix v1.0.2 fixes the local-package import failure seen in v1.0.1.

## What changed

The HeptadPair source package is a project-local package and is not guaranteed to be
installed into the active conda environment. v1.0.2 discovers `heptadpair/phase67.py`
from the supplied project root (supporting both `ROOT/heptadpair` and
`ROOT/src/heptadpair` layouts), prepends its parent to `PYTHONPATH`, performs an
explicit import preflight, and then runs the factorial analysis.

No feature definition, split, learner contract, target, or statistical analysis was
changed by this hotfix.

## Frozen inputs

- engineered: `results/features/heptadpair_mvp_features` -> 2214 features
- local register no-torsion: `results/features/two_center_topology_features_v064` -> 102
- ordered register no-torsion: `results/features/ordered_two_center_topology_features_v065` -> 118
- Phase67 locked splits: seeds 67-71
- reference replay: `results/phase67/final_frozen/model_summary.csv`

## Evaluated combinations

`EMPTY, E, O, L, EO, EL, OL, EOL`.

The scientific primary interpretation remains leave-one-block-out loss from the full
EOL model, supplemented by standalone, factorial interaction, and exact 3-block
Shapley summaries.

## Run

```bash
cd ~/BioData/hgq/Heptad/HeptadPair_v0_1
unzip -o HeptadPair_Phase81_FullFactorial_Ablation_v1_0_2_20260902.zip
chmod +x heptadpair_factorial_ablation_pack/*.sh
bash heptadpair_factorial_ablation_pack/run_phase81_factorial_ablation.sh "$PWD" \
  2>&1 | tee run_phase81_factorial_ablation_v1_0_2.log
```

A successful start must print lines like:

```text
[IMPORT PREFLIGHT PASS] heptadpair.phase67=.../heptadpair/phase67.py
[IMPORT PREFLIGHT PASS] KEYS=[...]
```

## Pack outputs

```bash
bash heptadpair_factorial_ablation_pack/pack_results.sh "$PWD"
```
