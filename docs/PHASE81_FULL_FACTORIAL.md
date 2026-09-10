# Phase81 full-factorial E/O/L ablation

Phase81 was added after the frozen primary CCNG1 benchmark to separate order-dependent incremental gains from block contribution.

Mean Spearman across the same five frozen outer splits:

| Subset | Dimensions | Mean Spearman |
|---|---:|---:|
| E | 2,214 | 0.603966 |
| O | 118 | 0.587648 |
| L | 102 | 0.526955 |
| E+O | 2,332 | 0.627549 |
| E+L | 2,316 | 0.623819 |
| O+L | 220 | 0.603125 |
| E+O+L | 2,434 | 0.632079 |

Primary leave-one-block-out Spearman losses from E+O+L:

- E: +0.028955
- O: +0.008261
- L: +0.004530

Supported interpretation: **conditional contribution in the final representation is E > O > L**. O is stronger than L as a stand-alone block, while the compact O+L representation (220 dimensions) nearly matches E (2,214 dimensions). The earlier E -> E+O -> E+O+L sequence is therefore an order-dependent incremental result, not an order-invariant ranking of intrinsic importance.

See `results/phase81_factorial_ablation_v1_0_2/` and `source_data/csv/Fig2_*` for the complete tables.
