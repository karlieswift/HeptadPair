# Phase-8.2C interpretation policy

1. The selected exact CCmax sequences have no direct 100%-identity experimental PDB pair hit from Phase-8.2B.
2. Any coordinates generated in Phase-8.2C are **modeled/predicted structural illustrations**, not experimental structures.
3. Case selection is retrospective visualization of the already frozen Section-3.4 retrieval result. It is not a new performance benchmark and must not be used to retune HeptadPair.
4. The primary main case is fixed before structure modeling:
   - intended pair: 4H3651-AA + 4H895-AA;
   - strongest sharing-member competitor: 4H3651-AA + 4H3365-AA;
   - frozen XGB score margin: +0.156695.
5. The same shared chain is placed as chain A in both modeling inputs so the target and competitor can be compared under a common visual convention.
6. The canonical heptad phase is inferred from sequence chemistry only to annotate a/d and e/g positions; any spatial salt bridge or hydrophobic contact claim must be verified on the coordinates.
7. AlphaFold/ColabFold confidence metrics are model-confidence diagnostics, not experimental evidence of binding or specificity.
