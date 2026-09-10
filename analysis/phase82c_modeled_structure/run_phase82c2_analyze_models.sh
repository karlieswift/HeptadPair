#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
PACK="$(cd "$(dirname "$0")" && pwd)"
python "$PACK/phase82c_analyze_models.py" "$ROOT" 2>&1 | tee "$ROOT/results/phase82c_modeled_structure/run_model_analysis.log"
