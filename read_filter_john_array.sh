#!/bin/bash
#SBATCH -J read_filter_john_array     # Job Name                # chmod +x /home/link/NB_EXP030/read_filter_john_array.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/NB_EXP030/read_filter_john_array.sh
#SBATCH --mem 8g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 00:60:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH --array=1-30
#SBATCH -o /g/steinmetz/link/logs/read_filter_john_array_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/read_filter_john_array_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail


set -euo pipefail

# Setup
input_dir="/g/steinmetz/project/to_lukas/NBEXP030"
output_dir="/g/steinmetz/link/Amplicon_barcode_analysis/NB_EXP030/john_read_filt"
thresholds="$output_dir/thresholds.tsv"
bam_list="$output_dir/bam_list.txt"

mkdir -p "$output_dir"

# Create the file with the bam paths as index
find "$input_dir" -type f -name '[IUL][L][1-4][1-4]_output_all_bc_umi.bam' | sort > "$bam_list"

# Get the BAM file for this task
bam_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$bam_list")
bam_base=$(basename "$bam_file")

# Parse replicate name
prefix=${bam_base%%_output_all_bc_umi.bam}
bin=${prefix:0:1}
sublib=${prefix:1:2}
sample=${prefix:3:1}
replicate_name="${bin}_${sublib}_${sample}"

# Lookup threshold
threshold=$(awk -v rep="$replicate_name" '$1 == rep {print $2}' "$thresholds")

# Skip if no threshold found
if [[ -z "$threshold" ]]; then
  echo "⚠️ No threshold for $replicate_name — skipping"
  exit 0
fi

# Output file paths
output_bam="$output_dir/${replicate_name}.bam"
counts_txt="$output_dir/${replicate_name}_umi_counts.txt"
counts_filt_txt="$output_dir/${replicate_name}_umi_counts_filt.txt"
results_tsv="$output_dir/${replicate_name}_results.tsv"

echo "---------------------------------------"
echo "Begin Processing: "
echo "BAM file:       $bam_file"
echo "Replicate name: $replicate_name"
echo "Threshold:      $threshold"
echo "Output BAM:     $output_bam"
echo "UMI counts:     $counts_txt"
echo "Filtered UMIs:  $counts_filt_txt"
echo "Results TSV:    $results_tsv"

# Load and use samtools
module purge
module load SAMtools

samtools view "$bam_file" | awk '
{
  cb = ""; umi = "";
  if (NF >= 18 && $18 ~ /^UB:Z:/ && $12 ~ /^CB:Z:/) {
    split($18, a, ":");
    split($12, b, ":");
    umi = a[3];
    cb = b[3];
  } else {
    for (i=NF; i>=12; i--) {
      if ($i ~ /^UB:Z:/) { split($i, a, ":"); umi = a[3]; }
      if ($i ~ /^CB:Z:/) { split($i, b, ":"); cb = b[3]; }
    }
  }
  if (length(umi) == 12 && cb != "") {
    print cb "|" umi;
  }
}' | sort | uniq -c > "$counts_txt"


echo "Finished generating $counts_txt"

awk -v T="$threshold" '$1 > T {print $2}' "$counts_txt" > "$counts_filt_txt"

echo "Finished generating $counts_filt_txt"
num_valid=$(wc -l < "$counts_filt_txt")
echo "UMIs passing threshold: $num_valid"

samtools view -h "$bam_file" | awk -v fname="$counts_filt_txt" '
BEGIN {
  while ((getline < fname) > 0) valid[$1]=1;
}
{
  if ($0 ~ /^@/) { print; next }

  cb = ""; umi = "";

  if (NF >= 18 && $18 ~ /^UB:Z:/ && $12 ~ /^CB:Z:/) {
    split($18, a, ":");
    split($12, b, ":");
    umi = a[3];
    cb = b[3];
  } else {
    for (i=NF; i>=12; i--) {
      if ($i ~ /^UB:Z:/) { split($i, a, ":"); umi = a[3]; }
      if ($i ~ /^CB:Z:/) { split($i, b, ":"); cb = b[3]; }
    }
  }

  key = cb "|" umi;
  if (length(umi) == 12 && (key in valid)) {
    print;
  }
}' | samtools view -b -o "$output_bam"


echo "Finished generating $output_bam"
samtools index "$output_bam"
echo "Finished indexing $output_bam"

# Load and run umi-tools
module purge
module load UMI-tools

umi_tools count \
  --per-gene \
  --gene-tag=GX \
  --umi-tag=UB \
  --extract-umi-method=tag \
  -I "$output_bam" \
  -S "$results_tsv"

echo "Finished generating $results_tsv"
echo "Finished processing $replicate_name"
echo "---------------------------------------"
