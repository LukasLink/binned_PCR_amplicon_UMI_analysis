#!/usr/bin/env bash
#SBATCH -J BCWQC_ARRAY_DCA                # chmod +x ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/bcwithqc_array_run.sh
#SBATCH -A lsteinme                       # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/bcwithqc_array_run.sh
#SBATCH --mem=40g
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH -t 96:00:00
#SBATCH --qos normal
#SBATCH --array=0-11
#SBATCH -o /g/steinmetz/link/logs/log_%x_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/log_%x_%A_%a.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=lukas.link@embl.de

set -euo pipefail

################################################################################
# USER OPTIONS
################################################################################

OUTPUT_FOLDER="/scratch/link/Amplicon_barcode_analysis/HepG2_dual_rep_DCA_bcwithqc"

BCWITHQC_CONFIG="/g/steinmetz/link/Amplicon_barcode_analysis/bcwithqc_test/bcwithqc_config_stagger.json"
STAR_INDEX="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_DCA/star_index/NOPE"
BCWITHQC_BIN="/home/link/.conda/envs/py313/bin/bcwithqc"

NUM_THREADS="${SLURM_CPUS_PER_TASK:-10}"
KEEP_INTERMEDIARY=false
VERBOSITY="-vv"

################################################################################
# END USER OPTIONS
################################################################################

timestamp() { date +"%T"; }
log() { echo "$(timestamp)   $*"; }

normalize_path() {
  local p="${1:-}"
  [[ -z "$p" ]] && { printf '%s' "$p"; return; }
  while [[ "$p" == *"//"* ]]; do p="${p//\/\//\/}"; done
  while [[ "$p" == */ && "$p" != "/" ]]; do p="${p%/}"; done
  printf '%s' "$p"
}

make_clean_dir() {
  mkdir -p "$1"
  normalize_path "$1"
}

check_file_exists() {
  local f="$1"
  local msg="${2:-Expected file not found}"
  if [[ ! -f "$f" ]]; then
    echo "ERROR: ${msg}: ${f}" >&2
    exit 1
  fi
}

check_dir_exists() {
  local d="$1"
  local msg="${2:-Expected directory not found}"
  if [[ ! -d "$d" ]]; then
    echo "ERROR: ${msg}: ${d}" >&2
    exit 1
  fi
}

OUTPUT_FOLDER="$(normalize_path "$OUTPUT_FOLDER")"
STAR_INDEX="$(normalize_path "$STAR_INDEX")"
BCWITHQC_CONFIG="$(normalize_path "$BCWITHQC_CONFIG")"

RUNS_DIR="$(make_clean_dir "$OUTPUT_FOLDER/bcwithqc_output")"
TMP_INPUT_DIR="$(normalize_path "$OUTPUT_FOLDER/tmp_inputs")"
LOGS_DIR="$(make_clean_dir "$OUTPUT_FOLDER/logs")"
MANIFEST="$(normalize_path "$OUTPUT_FOLDER/tmp_inputs_manifest.tsv")"

check_dir_exists "$TMP_INPUT_DIR" "Prepared symlink input directory not found; run bcwithqc_prepare_symlinks.sh first"
check_dir_exists "$STAR_INDEX" "STAR index folder not found"
check_file_exists "$BCWITHQC_CONFIG" "bcwithqc config not found"
check_file_exists "$MANIFEST" "Manifest not found; run bcwithqc_prepare_symlinks.sh first"

if [[ ! -x "$BCWITHQC_BIN" ]]; then
  echo "ERROR: bcwithqc executable is not executable or not found: $BCWITHQC_BIN" >&2
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" || -z "${SLURM_ARRAY_TASK_MIN:-}" || -z "${SLURM_ARRAY_TASK_MAX:-}" ]]; then
  echo "ERROR: This script must be submitted as a Slurm array job." >&2
  exit 1
fi

ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID"
ARRAY_TASK_MIN="$SLURM_ARRAY_TASK_MIN"
ARRAY_TASK_MAX="$SLURM_ARRAY_TASK_MAX"
ARRAY_STEP="${SLURM_ARRAY_TASK_STEP:-1}"

# This implementation assumes a contiguous array such as --array=0-2 or --array=1-3.
# For --array=0-2, zero_based_task_index is 0, 1, or 2 and total_array_tasks is 3.
if [[ "$ARRAY_STEP" -ne 1 ]]; then
  echo "ERROR: This script expects a contiguous Slurm array with step size 1, e.g. --array=0-2." >&2
  exit 1
fi

zero_based_task_index=$(( ARRAY_TASK_ID - ARRAY_TASK_MIN ))
total_array_tasks=$(( ARRAY_TASK_MAX - ARRAY_TASK_MIN + 1 ))

if (( zero_based_task_index < 0 || zero_based_task_index >= total_array_tasks )); then
  echo "ERROR: Invalid array task index calculation." >&2
  exit 1
fi

log "Using bcwithqc: $BCWITHQC_BIN"
log "Output folder: $OUTPUT_FOLDER"
log "Config: $BCWITHQC_CONFIG"
log "STAR index: $STAR_INDEX"
log "Threads per run: $NUM_THREADS"
log "Manifest: $MANIFEST"
log "Array task: $ARRAY_TASK_ID; zero-based index: $zero_based_task_index; total array tasks: $total_array_tasks"

################################################################################
# One bcwithqc + STAR + bcwithqc count run
################################################################################

run_bcwithqc_one_input_dir() {
  local input_dir="$1"
  local run_name="$2"

  local run_output_dir="$RUNS_DIR/$run_name"
  local star_dir="$run_output_dir/STAR_files"

  mkdir -p "$run_output_dir" "$star_dir"

  log "Starting run: $run_name"
  log "Run input dir: $input_dir"
  log "Run output dir: $run_output_dir"

  ##############################################################################
  # 1. bcwithqc preprocess
  ##############################################################################

  log "[$run_name] Starting bcwithqc preprocess"

  preprocess_cmd=(
    "$BCWITHQC_BIN" "preprocess"
    "$input_dir"
    "--config=$BCWITHQC_CONFIG"
    "--output-dir=$run_output_dir"
    "--threads=$NUM_THREADS"
  )

  if [[ -n "$VERBOSITY" ]]; then
    preprocess_cmd+=("$VERBOSITY")
  fi

  "${preprocess_cmd[@]}" 2>&1 | tee "$LOGS_DIR/${run_name}_01_preprocess.log"

  log "[$run_name] Finished bcwithqc preprocess"

  ##############################################################################
  # Find preprocessed FASTQ
  ##############################################################################

  mapfile -t preprocessed_fastqs < <(
    find "$run_output_dir" -maxdepth 1 \
      \( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \) \
      -type f | sort
  )

  if [[ "${#preprocessed_fastqs[@]}" -eq 0 ]]; then
    echo "ERROR: [$run_name] No preprocessed FASTQ found in: $run_output_dir" >&2
    exit 1
  fi

  if [[ "${#preprocessed_fastqs[@]}" -gt 1 ]]; then
    log "[$run_name] Multiple preprocessed FASTQs found:"
    printf '  %s\n' "${preprocessed_fastqs[@]}"
    echo "ERROR: Expected exactly one preprocessed FASTQ per run." >&2
    exit 1
  fi

  local preprocessed_fastq="${preprocessed_fastqs[0]}"

  log "[$run_name] Preprocessed FASTQ: $preprocessed_fastq"

  ##############################################################################
  # 2. STAR alignment
  ##############################################################################

  log "[$run_name] Starting STAR alignment"

  module purge
  ml STAR/2.7.11b-GCC-13.2.0

  local star_prefix="$star_dir/${run_name}_"

  star_cmd=(
    STAR
    --runThreadN "$NUM_THREADS"
    --genomeDir "$STAR_INDEX"
    --readFilesIn "$preprocessed_fastq"
    --outFileNamePrefix "$star_prefix"
    --outFilterMultimapNmax 1
    --outSAMtype BAM Unsorted
    --outSAMattributes NH HI AS nM GX GN
  )

  if [[ "$preprocessed_fastq" == *.gz ]]; then
    star_cmd+=(--readFilesCommand zcat)
  fi

  "${star_cmd[@]}" 2>&1 | tee "$LOGS_DIR/${run_name}_02_STAR.log"

  local aligned_bam="${star_prefix}Aligned.out.bam"
  check_file_exists "$aligned_bam" "[$run_name] Expected STAR BAM not found"

  log "[$run_name] Finished STAR alignment"
  log "[$run_name] STAR BAM: $aligned_bam"

  ##############################################################################
  # 3. bcwithqc count
  ##############################################################################

  log "[$run_name] Starting bcwithqc count"

  count_cmd=(
    "$BCWITHQC_BIN" "count"
    "$run_output_dir"
    "--STAR-output-dir=$star_dir"
    "--config=$BCWITHQC_CONFIG"
    "--output-dir=$run_output_dir"
    "--threads=$NUM_THREADS"
  )

  if [[ "$KEEP_INTERMEDIARY" == "true" ]]; then
    count_cmd+=("--keep-intermediary")
  fi

  if [[ -n "$VERBOSITY" ]]; then
    count_cmd+=("$VERBOSITY")
  fi

  "${count_cmd[@]}" 2>&1 | tee "$LOGS_DIR/${run_name}_03_count.log"

  log "[$run_name] Finished bcwithqc count"
  log "[$run_name] Completed run"
}

################################################################################
# Main execution: stride through prepared directories
################################################################################

# Read manifest, skipping header. Each entry is sample_name<TAB>input_dir<TAB>fastq.
declare -a sample_names=()
declare -a input_dirs=()

while IFS=$'\t' read -r sample_name input_dir fastq; do
  [[ "$sample_name" == "sample_name" ]] && continue
  [[ -z "$sample_name" ]] && continue
  sample_names+=("$sample_name")
  input_dirs+=("$input_dir")
done < "$MANIFEST"

total_samples="${#sample_names[@]}"

if [[ "$total_samples" -eq 0 ]]; then
  echo "ERROR: Manifest contains no samples: $MANIFEST" >&2
  exit 1
fi

log "Total prepared samples: $total_samples"

processed_any=false

for (( sample_index=zero_based_task_index; sample_index<total_samples; sample_index+=total_array_tasks )); do
  sample_name="${sample_names[$sample_index]}"
  input_dir="${input_dirs[$sample_index]}"

  check_dir_exists "$input_dir" "Prepared input directory missing for $sample_name"

  log "Array task $ARRAY_TASK_ID processing manifest index $sample_index: $sample_name"
  run_bcwithqc_one_input_dir "$input_dir" "$sample_name"
  processed_any=true
done

if [[ "$processed_any" == "false" ]]; then
  log "No samples assigned to array task $ARRAY_TASK_ID. Finishing successfully."
else
  log "Array task $ARRAY_TASK_ID completed all assigned samples."
fi
