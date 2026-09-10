# HeptadPair v8.9 GitHub release validation report

Release date: 2026-09-03

This compact GitHub/reviewer snapshot was rebuilt after the Phase81 full-factorial E/O/L revision and the Phase82 modeled structural-case revision were scientifically frozen. The release does not reopen model selection or retraining.

## Automated checks

- `scripts/repo_doctor.py --strict`: PASS; no required-file, cache, oversized-file or private-path issues.
- `scripts/verify_release.py`: PASS; all manifest hashes reproduce.
- `scripts/check_source_data.py`: PASS.
  - base-manuscript workbook: 16 sheets; the only formulas are the three frozen cumulative feature-accounting formulas in `Feature_Accounting!C2:C4`;
  - v8.9 revision workbook: 15 sheets; no formulas;
  - neither workbook contains known private/absolute filesystem identifiers.
- `tests/test_release_smoke.py`: 2/2 PASS.
- private-string scan across text and binary payloads: 0 hits for the development username, hostnames or absolute development paths.

## Reproduction checks

- `scripts/reproduce_revision_figures.py` reproduced the packaged PNGs for revised Figure 2, revised Figure 4 and Supplementary Figure S6 pixel-for-pixel in the validation environment.
- `scripts/compile_manuscripts.sh` compiled the repository LaTeX sources successfully:
  - main manuscript: 21 pages;
  - Supplementary Information: 12 pages.

## Release gate

The scientific/reviewer package is ready. `PUBLIC_RELEASE_GATE.json` intentionally keeps `public_release_ready=false` because license selection, final repository/archive identifiers, funding, detailed CRediT roles and final competing-interest wording remain author-controlled metadata.

No software/content license is inferred by this release.
