#!/bin/bash
#SBATCH -J read_filter_john     # Job Name                # chmod +x /home/link/NB_EXP030/read_filter_john.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/NB_EXP030/read_filter_john.sh
#SBATCH --mem 1g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 00:00:20             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/read_filter_john.out
#SBATCH -e /g/steinmetz/link/logs/read_filter_john.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

test_only=false

input_dir="/g/steinmetz/project/to_lukas/NBEXP030"
output_dir="/g/steinmetz/link/Amplicon_barcode_analysis/NB_EXP030/john_read_filt"
thresholds="$output_dir/thresholds.tsv"

mkdir -p "$output_dir"
find "$input_dir" -type f -name '[IUL][L][1-4][1-4]_output_all_bc_umi.bam' | sort > "$output_dir/bam_list.txt"


# Loop over matching BAM files
find "$input_dir" -type f -name '[IUL][L][1-4][1-4]_output_all_bc_umi.bam' | while read -r bam_file; do

  # Extract basename: e.g., LL32_output_all_bc_umi.bam
  bam_base=$(basename "$bam_file")

  prefix=${bam_base%%_output_all_bc_umi.bam}  # e.g., IL12
  # Extract parts of replicate name using pattern matching
  bin=${prefix:0:1}          # I or U or L
  sublib=${prefix:1:2}         # L1, L2, etc.
  sample=${prefix:3:1}          # final digit
  replicate_name="${bin}_${sublib}_${sample}"

  # Lookup threshold from TSV (column 2 where column 1 matches replicate_name)
  threshold=$(awk -v rep="$replicate_name" '$1 == rep {print $2}' "$thresholds")
  
  if [[ -z "$threshold" ]]; then
    echo "⚠️ WARNING: No threshold found for $replicate_name — skipping."
    continue
  fi

  # Define output file paths
  output_bam="$output_dir/${replicate_name}.bam"
  counts_txt="$output_dir/${replicate_name}_umi_counts.txt"
  counts_filt_txt="$output_dir/${replicate_name}_umi_counts_filt.txt"
  results_tsv="$output_dir/${replicate_name}_results.tsv"
  
  # Dry run output for validation
  echo "---------------------------------------"
  echo "Begin Processing: "
  echo "BAM file:       $bam_file"
  echo "Replicate name: $replicate_name"
  echo "Threshold:      $threshold"
  echo "Output BAM:     $output_bam"
  echo "UMI counts:     $counts_txt"
  echo "Filtered UMIs:  $counts_filt_txt"
  echo "Results TSV:    $results_tsv"
  
  module purge
  ml SAMtools
  samtools view "$bam_file" | awk '
  {
    umi = "";
    # Fast path
    if (NF >= 18 && $18 ~ /^UB:Z:/) {
      split($18, a, ":");
      umi = a[3];
    }
    # Slow path
    else {
      for (i=NF; i>=12; i--) {
        if ($i ~ /^UB:Z:/) {
          split($i, a, ":");
          umi = a[3];
          break;
        }
      }
    }
  
    if (length(umi) == 12) {
      print umi;
    }
  }' | sort | uniq -c > "$counts_txt"
  
  echo "Finished generating $counts_txt"
  
  awk -v T="$threshold" '$1 >= T {print $2}' "$counts_txt" > "$counts_filt_txt"
  
  echo "Finished generating $counts_filt_txt"
  num_valid=$(wc -l < "$counts_filt_txt")
  echo "UMIs passing threshold: $num_valid"
  
  samtools view -h "$bam_file" | \
  awk -v fname="$counts_filt_txt" '
  BEGIN {
    while ((getline < fname) > 0) valid[$1]=1;
  }
  {
    if ($0 ~ /^@/) { print; next }
  
    # Fast path: try field 18 directly
    if (NF >= 18 && $18 ~ /^UB:Z:/) {
      split($18, a, ":");
      if (a[3] in valid) {
        print;
      }
    }
    else {
      # Slow path: check all fields in reverse
      for (i=NF; i>=12; i--) {
        if ($i ~ /^UB:Z:/) {
          split($i, a, ":");
          if (a[3] in valid) {
            print;
          }
          break;  # stop after first UB found
        }
      }
    }
  }' | samtools view -b -o "$output_bam"
  echo "Finished generating $output_bam"
  samtools index "$output_bam"
  echo "Finished indexing $output_bam"
  
  module purge
  ml UMI-tools
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

done








# bam_file="/g/steinmetz/project/to_lukas/NBEXP030/LL32_output_all/LL32_output_all_bc_umi.bam"
# output_dir="/g/steinmetz/link/Amplicon_barcode_analysis/NB_EXP030/test"
# output_bam="$output_dir/out.bam"
# counts_txt="$output_dir/umi_counts.txt"
# counts_filt_txt="$output_dir/umi_counts_filt.txt"
# threshold=1000 #set the threshold
# results_tsv="$output_dir/results.tsv"
# 
# mkdir -p "$output_dir"
# 
# 
# 
# module purge
# ml SAMtools
# samtools view "$bam_file" | awk '
# {
#   umi = "";
#   # Fast path
#   if (NF >= 18 && $18 ~ /^UB:Z:/) {
#     split($18, a, ":");
#     umi = a[3];
#   }
#   # Slow path
#   else {
#     for (i=NF; i>=12; i--) {
#       if ($i ~ /^UB:Z:/) {
#         split($i, a, ":");
#         umi = a[3];
#         break;
#       }
#     }
#   }
# 
#   if (length(umi) == 12) {
#     print umi;
#   }
# }' | sort | uniq -c > "$counts_txt"
# 
# awk -v T="$threshold" '$1 >= T {print $2}' "$counts_txt" > "$counts_filt_txt"
# 
# samtools view -h "$bam_file" | \
# awk -v fname="$counts_filt_txt" '
# BEGIN {
#   while ((getline < fname) > 0) valid[$1]=1;
# }
# {
#   if ($0 ~ /^@/) { print; next }
# 
#   # Fast path: try field 18 directly
#   if (NF >= 18 && $18 ~ /^UB:Z:/) {
#     split($18, a, ":");
#     if (a[3] in valid) {
#       print;
#     }
#   }
#   else {
#     # Slow path: check all fields in reverse
#     for (i=NF; i>=12; i--) {
#       if ($i ~ /^UB:Z:/) {
#         split($i, a, ":");
#         if (a[3] in valid) {
#           print;
#         }
#         break;  # stop after first UB found
#       }
#     }
#   }
# }' | samtools view -b -o "$output_bam"
# 
# samtools index "$output_bam"
# 
# module purge
# ml UMI-tools
# umi_tools count \
#   --per-gene \
#   --gene-tag=GX \
#   --umi-tag=UB \
#   --extract-umi-method=tag \
#   -I "$output_bam" \
#   -S "$results_tsv"

