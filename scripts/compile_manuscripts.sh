#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
cd "$ROOT/manuscript"
pdflatex -interaction=nonstopmode -halt-on-error HeptadPair_Main_Manuscript_v8.9_PHASE81_PHASE82_20260903.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error HeptadPair_Main_Manuscript_v8.9_PHASE81_PHASE82_20260903.tex >/dev/null
cd "$ROOT/supplementary"
pdflatex -interaction=nonstopmode -halt-on-error HeptadPair_Supplementary_Information_v8.9_PHASE81_PHASE82_20260903.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error HeptadPair_Supplementary_Information_v8.9_PHASE81_PHASE82_20260903.tex >/dev/null
echo "[PASS] manuscript and SI compiled"
