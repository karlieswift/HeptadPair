# HeptadPair Phase82D - matched-orientation off-target control

Purpose: add a **secondary** same-orientation control for the fixed Phase82C case.
It does not replace the strongest overall competitor and does not retrain/retune
HeptadPair.

Selection is locked before structure modeling:
1. Keep target `4H3651-AA + 4H895-AA` fixed from Phase82C.
2. Use all heterodimer off-targets in the same official set that share `4H3651-AA`.
3. Order them by the existing frozen `xgb_x2_1250_final_shared_2434` score.
4. Model all remaining candidates with the exact Phase82C AF2-Multimer protocol.
5. Choose the highest-scoring candidate whose 25/25 modeled orientations are parallel.
6. If no such candidate exists, report that result; do not loosen the rule post hoc.

Run Phase82D0 under `bioh200`, Phase82D1 under `colabfold82c`, and Phase82D2 under `bioh200`.
