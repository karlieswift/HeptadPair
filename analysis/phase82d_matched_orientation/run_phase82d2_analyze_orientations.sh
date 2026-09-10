#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
PACK="$(cd "$(dirname "$0")" && pwd)"
python "$PACK/phase82d_analyze_orientations.py" "$ROOT" 2>&1 | tee "$ROOT/results/phase82d_matched_orientation/run_orientation_analysis.log"
