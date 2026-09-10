#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
PYTHON="${PYTHON:-python}"
OUT="$ROOT/results/phase82_rcsb_screen"
rm -rf "$OUT"
"$PYTHON" "$(dirname "$0")/selftest_phase82b_payload.py"
"$PYTHON" "$(dirname "$0")/phase82b_rcsb_exact_sequence_screen.py" "$ROOT" "${@:2}"
