#!/bin/bash
#SBATCH -J split_input_2     # Job Name                # chmod +x /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/split_input_for_fake_untreated.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/split_input_for_fake_untreated.sh
#SBATCH --mem 36g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 04:00:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/log_split_input_2.out
#SBATCH -e /g/steinmetz/link/logs//log_split_input_2.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################
# Input and output folder paths
INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/Dual_rep_quart_rush"
OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/fake_untreated_2"
percentages=("75" "50" "25")

################################################################################
# END OF USER OPTIONS
################################################################################

ml seqtk
ml pigz   # if available

# Create QC_filtered_split directory (or reuse QC_filtered_subsampled if you like)
QC_DIR="$OUTPUT_FOLDER/QC_filtered_split"
mkdir -p "$QC_DIR"
echo "Split directory created successfully!"

# Split .fastq.gz files starting with I_ in QC_filtered directory
if [ -d "$INPUT_FOLDER/QC_filtered" ]; then
    for file in "$INPUT_FOLDER/QC_filtered"/I_*.fastq.gz; do
        # Handle case where glob matches nothing
        [ -e "$file" ] || continue

        fname=$(basename "$file")

        # Output file names
        outI="$QC_DIR/$fname"
        outL="$QC_DIR/${fname/I_/L_}"
        outU="$QC_DIR/${fname/I_/U_}"

        echo "Splitting $fname into:"
        echo "  $outI"
        echo "  $outL"
        echo "  $outU"

        zcat "$file" | awk -v oI="$outI" -v oL="$outL" -v oU="$outU" -v fname="$fname" '
            BEGIN {
                srand(444);  # fixed seed for reproducibility (change/remove if you want)
                cmdI = "pigz -p 4 > " oI;  # adjust -p 4 to your core count
                cmdL = "pigz -p 4 > " oL;
                cmdU = "pigz -p 4 > " oU;
                readN = 0;
            }
            {
                # Every 4 lines = one FASTQ read
                if (NR % 4 == 1) {
                    readN++;
                    r = rand();

                    # Progress messages
                    if (readN == 1000 || readN == 10000 || readN == 100000 || readN % 1000000 == 0) {
                        # Print to stderr so it doesn’t interfere with pipelines
                        printf("File %s: processed %d reads\n", fname, readN) > "/dev/stderr";
                    }
                }

                if (r < 1.0/3.0) {
                    print | cmdI;
                } else if (r < 2.0/3.0) {
                    print | cmdL;
                } else {
                    print | cmdU;
                }
            }
            END {
                close(cmdI);
                close(cmdL);
                close(cmdU);
            }
        '
    done
fi

echo "All I_* fastq files split successfully!"