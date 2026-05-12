#!/usr/bin/env bash 
#SBATCH -J BCWQC                 # chmod +x ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/bcwithqc_pipeline.sh
#SBATCH -A lsteinme              # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/bcwithqc_pipeline.sh
#SBATCH --mem 128g               #  --dependency=afterok:46151850
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH --qos normal
#SBATCH -o /g/steinmetz/link/logs/log_%x_%A.out
#SBATCH -e /g/steinmetz/link/logs/log_%x_%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lukas.link@embl.de

set -euo pipefail

################################################################################
# USER OPTIONS
################################################################################

INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/1/HepG2_dual_rep_PA_subsample/subsample_10/QC_filtered/bcwithqc_test"
OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/bcwithqc_test/stagger_test"

BCWITHQC_CONFIG="/g/steinmetz/link/Amplicon_barcode_analysis/bcwithqc_test/bcwithqc_config_stagger.json"

# Existing STAR index
STAR_INDEX="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/1/HepG2_dual_rep_PA_subsample/subsample_10/star_index/NOPE"

# Use the working bcwithqc executable directly.
# Later this can become simply: bcwithqc
BCWITHQC_BIN="/home/link/.conda/envs/py313/bin/bcwithqc"

NUM_THREADS=10

# If true, pass the complete INPUT_FOLDER to one bcwithqc run.
# If false, create one lightweight symlink-input directory per FASTQ file.
ENTIRE_FOLDER_AS_ONE=false

# File patterns to process when ENTIRE_FOLDER_AS_ONE=false
FASTQ_PATTERNS=("*.fastq.gz" "*.fq.gz" "*.fastq" "*.fq")

# Keep bcwithqc intermediary files
KEEP_INTERMEDIARY=true

# Additional bcwithqc verbosity option.
# Example: "-v", "-vv", or "".
VERBOSITY="-vv"

################################################################################
# END USER OPTIONS
################################################################################


timestamp() {
  date +"%T"
}

log() {
  echo "$(timestamp)   $*"
}

normalize_path() {
  local p="${1:-}"

  [[ -z "$p" ]] && { printf '%s' "$p"; return; }

  while [[ "$p" == *"//"* ]]; do
    p="${p//\/\//\/}"
  done

  while [[ "$p" == */ && "$p" != "/" ]]; do
    p="${p%/}"
  done

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

INPUT_FOLDER="$(normalize_path "$INPUT_FOLDER")"
OUTPUT_FOLDER="$(normalize_path "$OUTPUT_FOLDER")"
STAR_INDEX="$(normalize_path "$STAR_INDEX")"
BCWITHQC_CONFIG="$(normalize_path "$BCWITHQC_CONFIG")"

check_dir_exists "$INPUT_FOLDER" "Input folder not found"
check_dir_exists "$STAR_INDEX" "STAR index folder not found"
check_file_exists "$BCWITHQC_CONFIG" "bcwithqc config not found"

if [[ ! -x "$BCWITHQC_BIN" ]]; then
  echo "ERROR: bcwithqc executable is not executable or not found: $BCWITHQC_BIN" >&2
  exit 1
fi

RUNS_DIR="$(make_clean_dir "$OUTPUT_FOLDER/bcwithqc_output")"
TMP_INPUT_DIR="$(make_clean_dir "$OUTPUT_FOLDER/tmp_inputs")"
LOGS_DIR="$(make_clean_dir "$OUTPUT_FOLDER/logs")"

log "Using bcwithqc: $BCWITHQC_BIN"
log "Input folder: $INPUT_FOLDER"
log "Output folder: $OUTPUT_FOLDER"
log "Config: $BCWITHQC_CONFIG"
log "STAR index: $STAR_INDEX"
log "Threads: $NUM_THREADS"


################################################################################
# Optional: rename .txt.gz files to .fastq.gz
################################################################################

for f in "$INPUT_FOLDER"/*.txt.gz; do
  [[ -e "$f" ]] || continue
  new_name="${f%.txt.gz}.fastq.gz"
  log "Renaming: $f -> $new_name"
  mv -- "$f" "$new_name"
done


################################################################################
# Collect input files
################################################################################

collect_fastq_files() {
  local files=()

  for pattern in "${FASTQ_PATTERNS[@]}"; do
    while IFS= read -r -d '' f; do
      files+=("$f")
    done < <(find "$INPUT_FOLDER" -maxdepth 1 -type f -name "$pattern" -print0)
  done

  printf '%s\n' "${files[@]}"
}


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
  # bcwithqc typically writes a sans_bc*.fq / *.fastq file.
  # We search robustly instead of hardcoding one exact filename.

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
    echo "ERROR: For this first script, expected exactly one preprocessed FASTQ per run." >&2
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
# Main execution
################################################################################

if [[ "$ENTIRE_FOLDER_AS_ONE" == "true" ]]; then
  log "ENTIRE_FOLDER_AS_ONE=true"
  log "Passing complete input folder to one bcwithqc run"

  run_bcwithqc_one_input_dir "$INPUT_FOLDER" "all_files"

else
  log "ENTIRE_FOLDER_AS_ONE=false"
  log "Creating one symlink-based input directory per FASTQ"

  mapfile -t fastq_files < <(collect_fastq_files)

  if [[ "${#fastq_files[@]}" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in $INPUT_FOLDER" >&2
    exit 1
  fi

  log "Found ${#fastq_files[@]} FASTQ files"

  for fastq in "${fastq_files[@]}"; do
    filename="$(basename "$fastq")"

    sample_name="$filename"
    sample_name="${sample_name%.fastq.gz}"
    sample_name="${sample_name%.fq.gz}"
    sample_name="${sample_name%.fastq}"
    sample_name="${sample_name%.fq}"

    sample_input_dir="$TMP_INPUT_DIR/$sample_name"
    mkdir -p "$sample_input_dir"

    # Clean previous symlink if rerunning.
    rm -f "$sample_input_dir/$filename"

    # Lightweight fake directory: contains only a symlink to the original FASTQ.
    ln -s "$fastq" "$sample_input_dir/$filename"

    log "Prepared symlink input dir for $sample_name:"
    log "  $sample_input_dir/$filename -> $fastq"

    run_bcwithqc_one_input_dir "$sample_input_dir" "$sample_name"
  done
fi

log "Processing complete."