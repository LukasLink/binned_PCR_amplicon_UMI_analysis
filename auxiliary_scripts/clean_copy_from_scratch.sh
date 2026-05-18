#!/usr/bin/env bash
#SBATCH -J clean_copy               # chmod +x ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/clean_copy.sh
#SBATCH -A lsteinme                      # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/clean_copy.sh
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

SRC="/scratch/link/Amplicon_barcode_analysis"
DST="/g/steinmetz/link/Amplicon_barcode_analysis"

case "$DST/" in "$SRC"/*)
  echo "ERROR: destination must not be inside source directory."
  exit 1
  ;;
esac

move_unique() {
  local src="$1" dst="$2" dir base stem ext candidate n
  dir=$(dirname "$dst")
  base=$(basename "$dst")
  stem="${base%.*}"
  ext=".${base##*.}"
  [[ "$stem" == "$base" ]] && ext=""
  candidate="$dst"
  n=1

  while [[ -e "$candidate" ]]; do
    candidate="$dir/${stem}_${n}${ext}"
    ((n++))
  done

  mv -- "$src" "$candidate"
}

mapfile -d '' INTER_DIRS < <(find "$SRC" -depth -type d -name "intermediary_files" -print0)

for idir in "${INTER_DIRS[@]}"; do
  [[ -d "$idir" ]] || continue

  parent=$(dirname "$idir")
  logdir="$parent/logs"
  mkdir -p "$logdir"

  while IFS= read -r -d '' f; do
    name=$(basename "$f")

    if [[ "$name" == *.out && "$name" != STAR_* ]]; then
      name="STAR_$name"
    fi

    move_unique "$f" "$logdir/$name"
  done < <(find "$idir" -type f \( -name "*.log" -o -name "*.out" \) -print0)

  rm -rf -- "$idir"
done

mkdir -p "$DST"
rsync -a -- "$SRC"/ "$DST"/