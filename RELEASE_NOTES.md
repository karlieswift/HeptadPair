# Release notes — v8.9 Phase81 + Phase82

This release supersedes the earlier manuscript-stage GitHub snapshot for the purposes of the current revision.

## Added

- complete Phase81 full-factorial E/O/L analysis and replay audits;
- revised Figure 2 source tables and reconstruction code;
- Phase82 score-first structural-case selection, exact experimental-PDB screen, AlphaFold2-Multimer model ensemble outputs and matched-orientation sensitivity controls;
- revised Figure 4 and Supplementary Figure S6 structural panels;
- public-safe source-data workbook with repository-relative paths;
- portable `repo_doctor.py`, `verify_release.py`, source-data checks and revised-figure reconstruction;
- CI workflow for release verification.

## Interpretation changes

- E/O/L importance is no longer inferred from the sequential E -> E+O -> E+O+L path. The supported conditional ranking is E > O > L.
- The structural example is explicitly described as modeled/predicted, because the exact experimental-PDB screen found no common direct structure.
- The intended partner is not claimed to win because it simply folds better: hard and matched-orientation off-targets are structurally plausible and can have denser modeled contacts.
