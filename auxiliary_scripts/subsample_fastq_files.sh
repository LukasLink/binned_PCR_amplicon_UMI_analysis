#!/bin/bash
#SBATCH -J subsample_fastq_files     # Job Name                # chmod +x /g/steinmetz/battisti/Data_analyses/ampliconseq/subsample_fastq_files.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /g/steinmetz/battisti/Data_analyses/ampliconseq/subsample_fastq_files.sh
#SBATCH --mem 128g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 12                   # Number of CPUs
#SBATCH -t 24:05:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/battisti/Data_analyses/ampliconseq/log_subsample_fastq_files.out
#SBATCH -e /g/steinmetz/battisti/Data_analyses/ampliconseq/log_subsample_fastq_files.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################
# Input and output folder paths
INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/RAW/Dual_rep_quart_rush"
OUTPUT_FOLDER="/g/steinmetz/battisti/Data_analyses/ampliconseq/Dual_rep_quart_rush"
percentages=("75" "50" "25")

################################################################################
# END OF USER OPTIONS
################################################################################

# Extract the input directory name from the path
INPUT_DIR_NAME=$(basename "$INPUT_FOLDER")

# Create the main output directory with the suffix "_subsample"
NEW_DIR="$OUTPUT_FOLDER/${INPUT_DIR_NAME}_subsample"
mkdir -p "$NEW_DIR"  # Create the directory, -p ensures no error if it already exists

# List of percentages to create subdirectories for


# Loop through the percentages and create subdirectories
for percent in "${percentages[@]}"; do
    SUBDIR="$NEW_DIR/subsample_${percent}"
    mkdir -p "$SUBDIR"  # Create the subdirectory for the percentage

    # Copy the ref directory to each subdirectory if it exists
    if [ -d "$INPUT_FOLDER/ref" ]; then
        cp -r "$INPUT_FOLDER/ref" "$SUBDIR/"
    fi

    # Copy files from the rds directory (excluding files with "MAUDE" in the name) to each subdirectory
    if [ -d "$INPUT_FOLDER/rds" ]; then
        # Find files in rds that don't contain "MAUDE" in the name and copy them to the subdirectory
        find "$INPUT_FOLDER/rds" -type f ! -name "*MAUDE*" -exec cp {} "$SUBDIR/" \;
    fi
    # Create QC_filtered_subsampled directory
    QC_DIR="$SUBDIR/QC_filtered_subsampled"
    mkdir -p "$QC_DIR"
    echo "Subdirectories created successfully!"
    
    # Subsample .fastq.gz or .txt.gz files in QC_filtered directory
    if [ -d "$INPUT_FOLDER/QC_filtered" ]; then
        for file in "$INPUT_FOLDER/QC_filtered"/*; do
            # Check for .fastq.gz files and subsample them
            if [[ "$file" == *.fastq.gz ]]; then
                # Calculate the number of reads to subsample
                TOTAL_READS=$(zcat "$file" | wc -l)
                READS_TO_SAMPLE=$((TOTAL_READS * percent / 100))

                # Use seqtk to randomly subsample the file and save it in the subsampled directory
                seqtk sample "$file" "$READS_TO_SAMPLE" | gzip > "$QC_DIR/$(basename "$file")"

            # Check for .txt.gz files and subsample them
            elif [[ "$file" == *.txt.gz ]]; then
                # Calculate the number of lines to subsample
                TOTAL_LINES=$(zcat "$file" | wc -l)
                LINES_TO_SAMPLE=$((TOTAL_LINES * percent / 100))

                # Use shuf to randomly subsample the lines and save it in the subsampled directory
                zcat "$file" | shuf -n "$LINES_TO_SAMPLE" | gzip > "$QC_DIR/$(basename "$file")"
            fi
        done
    fi
    echo "subsampled files created successfully!"
done
echo "All Done!"

