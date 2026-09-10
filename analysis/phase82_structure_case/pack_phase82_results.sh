#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="HeptadPair_Phase82_StructureCase_RESULTS_v1_2_0_${STAMP}.tar.gz"
cd "$ROOT"
items=()
[[ -d results/phase82_structure_case_preflight ]] && items+=(results/phase82_structure_case_preflight)
[[ -d results/phase82_rcsb_screen ]] && items+=(results/phase82_rcsb_screen)
[[ -f run_phase82a_structure_case_preflight_v1_2_0.log ]] && items+=(run_phase82a_structure_case_preflight_v1_2_0.log)
[[ -f run_phase82b_rcsb_exact_sequence_screen_v1_2_0.log ]] && items+=(run_phase82b_rcsb_exact_sequence_screen_v1_2_0.log)
if [[ ${#items[@]} -eq 0 ]]; then echo "[ERROR] no Phase82 outputs found" >&2; exit 2; fi
tar -czf "$OUT" "${items[@]}"
sha256sum "$OUT" > "$OUT.sha256"
echo "[OK] $ROOT/$OUT"
echo "[OK] $ROOT/$OUT.sha256"
