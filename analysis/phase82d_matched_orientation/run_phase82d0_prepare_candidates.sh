#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
PACK="$(cd "$(dirname "$0")" && pwd)"
python "$PACK/phase82d_prepare_candidates.py" "$ROOT"
