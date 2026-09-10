#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$PWD}"
ROOT="$(cd "$ROOT" && pwd)"
OUT="$ROOT/results/phase82c_modeled_structure"
INPUT="$OUT/inputs/main_pair_batch.fasta"
PRED="$OUT/colabfold_main"
GPU="${PHASE82C_GPU:-0}"
SEEDS="${PHASE82C_NUM_SEEDS:-5}"
MODELS="${PHASE82C_NUM_MODELS:-5}"
RECYCLES="${PHASE82C_NUM_RECYCLES:-20}"
RSEED="${PHASE82C_RANDOM_SEED:-20260902}"
if [[ ! -s "$INPUT" ]]; then echo "Missing $INPUT; run Phase82C0 first" >&2; exit 2; fi
if ! command -v colabfold_batch >/dev/null 2>&1; then
  echo "[BLOCK] colabfold_batch not found in PATH." >&2
  echo "Activate a ColabFold environment and rerun. Inputs are already frozen at $INPUT" >&2
  exit 3
fi
mkdir -p "$PRED"
HELP="$(colabfold_batch --help 2>&1 || true)"
args=("$INPUT" "$PRED")
if grep -q -- '--msa-mode' <<<"$HELP"; then args+=(--msa-mode single_sequence); fi
if grep -q -- 'alphafold2_multimer_v3' <<<"$HELP"; then args+=(--model-type alphafold2_multimer_v3); fi
if grep -q -- '--pair-mode' <<<"$HELP"; then args+=(--pair-mode unpaired); fi
if grep -q -- '--num-recycle' <<<"$HELP"; then args+=(--num-recycle "$RECYCLES"); fi
if grep -q -- '--num-seeds' <<<"$HELP"; then args+=(--num-seeds "$SEEDS"); fi
if grep -q -- '--random-seed' <<<"$HELP"; then args+=(--random-seed "$RSEED"); fi
if grep -q -- '--num-models' <<<"$HELP"; then args+=(--num-models "$MODELS"); fi
if grep -q -- '--num-relax' <<<"$HELP"; then args+=(--num-relax 0); fi
if grep -q -- '--rank' <<<"$HELP" && grep -q 'multimer' <<<"$HELP"; then args+=(--rank multimer); fi
if grep -q -- '--overwrite-existing-results' <<<"$HELP"; then args+=(--overwrite-existing-results); fi
printf '[RUN] CUDA_VISIBLE_DEVICES=%s colabfold_batch' "$GPU"
printf ' %q' "${args[@]}"
printf '\n'
CUDA_VISIBLE_DEVICES="$GPU" colabfold_batch "${args[@]}" 2>&1 | tee "$OUT/run_colabfold_main.log"
echo "[OK] $PRED"
