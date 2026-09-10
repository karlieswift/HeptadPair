#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
PYTHON="${PYTHON:-python}"
OUT="$ROOT/results/phase82_structure_case_preflight"
rm -rf "$OUT"
"$PYTHON" "$(dirname "$0")/selftest_phase82a.py"
"$PYTHON" "$(dirname "$0")/phase82a_structure_case_preflight.py" "$ROOT" "${@:2}"
