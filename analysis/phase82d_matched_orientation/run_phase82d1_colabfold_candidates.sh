#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
IN="$ROOT/results/phase82d_matched_orientation/matched_parallel_candidate_batch.fasta"
OUT="$ROOT/results/phase82d_matched_orientation/colabfold_candidates"
GPU="${PHASE82D_GPU:-0}"
NSEEDS="${PHASE82D_NUM_SEEDS:-5}"
NMODELS="${PHASE82D_NUM_MODELS:-5}"
[[ -s "$IN" ]] || { echo "Missing/empty $IN; run Phase82D0 first" >&2; exit 2; }
command -v colabfold_batch >/dev/null 2>&1 || { echo "colabfold_batch not found in active environment" >&2; exit 2; }
mkdir -p "$OUT"
echo "[RUN] CUDA_VISIBLE_DEVICES=$GPU colabfold_batch $IN $OUT"
CUDA_VISIBLE_DEVICES="$GPU" colabfold_batch "$IN" "$OUT" \
  --msa-mode single_sequence \
  --model-type alphafold2_multimer_v3 \
  --pair-mode unpaired \
  --num-recycle 20 \
  --num-seeds "$NSEEDS" \
  --random-seed 20260902 \
  --num-models "$NMODELS" \
  --num-relax 0 \
  --rank multimer \
  --overwrite-existing-results \
  2>&1 | tee "$ROOT/results/phase82d_matched_orientation/run_colabfold_candidates.log"
echo "[OK] $OUT"
