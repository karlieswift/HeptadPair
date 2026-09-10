# Phase82 partner-specific modeled structural illustration

The intended target and off-target were chosen from frozen final E+O+L XGB retrieval scores **before** structure modeling.

| Case | Frozen score | Orientation ensemble | Mean ipTM | Median contacts | Median a/d | Median e/g opp/same |
|---|---:|---|---:|---:|---:|---:|
| 4H3651 + 4H895 intended target | -0.510627 | 25/25 parallel | 0.7928 | 36 | 18 | 6 / 0 |
| 4H3651 + 4H3365 strongest off-target | -0.667322 | 25/25 antiparallel | 0.7808 | 44 | 21 | 5 / 1 |
| 4H3651 + 4H140 strict parallel control | -1.576181 | 25/25 parallel | 0.8184 | 46 | 22 | 5 / 3 |

The strongest off-target is not a trivial non-folding negative, and the strict parallel control has higher modeled complex confidence plus denser overall/a-d contact packing than the intended target yet receives a much lower HeptadPair score. The supported interpretation is partner-specific discrimination beyond generic modeled foldability, interface size or simple a/d-core packing.

No exact experimental PDB structure was found for the selected target/competitor. ipTM/pTM are model-confidence values rather than affinity measurements; e/g patterns are mechanistic illustrations rather than causal proof.
