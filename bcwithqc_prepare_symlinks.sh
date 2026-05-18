#!/usr/bin/env bash
#SBATCH -J BCWQC_PREP                # chmod +x ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/bcwithqc_prepare_symlinks.sh
#SBATCH -A lsteinme                      # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/bcwithqc_prepare_symlinks.sh
#SBATCH --mem=1g
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -t 00:01:00
#SBATCH --qos normal
#SBATCH -o /g/steinmetz/link/logs/log_%x_%j.out
#SBATCH -e /g/steinmetz/link/logs/log_%x_%j.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=lukas.link@embl.de

set -euo pipefail

################################################################################
# USER OPTIONS
################################################################################

INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_DCA/QC_filtered"
OUTPUT_FOLDER="/scratch/link/Amplicon_barcode_analysis/HepG2_dual_rep_DCA_bcwithqc"

FASTQ_PATTERNS=("*.fastq.gz" "*.fq.gz" "*.txt.gz" "*.fastq" "*.fq" "*.txt")

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

check_dir_exists() {
  local d="$1"
  local msg="${2:-Expected directory not found}"
  if [[ ! -d "$d" ]]; then
    echo "ERROR: ${msg}: ${d}" >&2
    exit 1
  fi
}

make_clean_dir() {
  mkdir -p "$1"
  normalize_path "$1"
}

collect_fastq_files() {
  local files=()
  for pattern in "${FASTQ_PATTERNS[@]}"; do
    while IFS= read -r -d '' f; do
      files+=("$f")
    done < <(find "$INPUT_FOLDER" -maxdepth 1 -type f -name "$pattern" -print0)
  done
  printf '%s\n' "${files[@]}" | sort -V
}

INPUT_FOLDER="$(normalize_path "$INPUT_FOLDER")"
OUTPUT_FOLDER="$(normalize_path "$OUTPUT_FOLDER")"
TMP_INPUT_DIR="$(make_clean_dir "$OUTPUT_FOLDER/tmp_inputs")"
MANIFEST="$(normalize_path "$OUTPUT_FOLDER/tmp_inputs_manifest.tsv")"

check_dir_exists "$INPUT_FOLDER" "Input folder not found"

log "Input folder: $INPUT_FOLDER"
log "Output folder: $OUTPUT_FOLDER"
log "Temporary symlink input folder: $TMP_INPUT_DIR"
log "Manifest: $MANIFEST"

################################################################################
# Create one symlink input directory per FASTQ
################################################################################

mapfile -t fastq_files < <(collect_fastq_files)

if [[ "${#fastq_files[@]}" -eq 0 ]]; then
  echo "ERROR: No FASTQ files found in $INPUT_FOLDER" >&2
  exit 1
fi

log "Found ${#fastq_files[@]} FASTQ files"

# Recreate the manifest on every preparation run.
printf 'sample_name\tinput_dir\tfastq\n' > "$MANIFEST"

for fastq in "${fastq_files[@]}"; do
  filename="$(basename "$fastq")"

  sample_name="$filename"
  sample_name="${sample_name%.fastq.gz}"
  sample_name="${sample_name%.fq.gz}"
  sample_name="${sample_name%.fastq}"
  sample_name="${sample_name%.fq}"

  sample_input_dir="$TMP_INPUT_DIR/$sample_name"
  mkdir -p "$sample_input_dir"

  # Keep directory clean in case the script is rerun.
  find "$sample_input_dir" -mindepth 1 -maxdepth 1 -type l -delete
  rm -f "$sample_input_dir/$filename"

  ln -s "$fastq" "$sample_input_dir/$filename"

  printf '%s\t%s\t%s\n' "$sample_name" "$sample_input_dir" "$fastq" >> "$MANIFEST"

  log "Prepared $sample_name"
  log "  $sample_input_dir/$filename -> $fastq"
done

log "Preparation complete. Created ${#fastq_files[@]} symlink input directories."
