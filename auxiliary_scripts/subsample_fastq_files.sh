#!/bin/bash
#SBATCH -J subsample_fastq_files_1     # Job Name                # chmod +x /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/subsample_fastq_files.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/subsample_fastq_files.sh
#SBATCH --mem 36g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 07:00:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/log_subsample_fastq_files_1.out
#SBATCH -e /g/steinmetz/link/logs//log_subsample_fastq_files_1.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################
# Input and output folder paths
# path to the master folder for this analysis example HepG2_dual_rep_PA
INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_PA"
# Path to where the ouput folder should be created, example Amplicon_barcode_analysis
OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/1"
percentages=("90" "80" "70" "60" "50" "40" "30" "20" "10")

################################################################################
# END OF USER OPTIONS
################################################################################

ml seqtk

# Extract the input directory name from the path
INPUT_DIR_NAME=$(basename "$INPUT_FOLDER")

# Create the main output directory with the suffix "_subsample"
NEW_DIR="$OUTPUT_FOLDER/${INPUT_DIR_NAME}_subsample"
mkdir -p "$NEW_DIR"  # Create the directory, -p ensures no error if it already exists

# List of percentages to create subdirectories for


# Loop through the percentages and create subdirectories
for percent in "${percentages[@]}"; do
    # Ensure percent is treated as a number
    percent=$((percent))  # This forces the string to be treated as a number
    percent_decimal=$(echo "$percent / 100" | bc -l)
    
    SUBDIR="$NEW_DIR/subsample_${percent}"
    mkdir -p "$SUBDIR"  # Create the subdirectory for the percentage

    # Copy the ref directory to each subdirectory if it exists
    if [ -d "$INPUT_FOLDER/ref" ]; then
        cp -r "$INPUT_FOLDER/ref" "$SUBDIR/"
    fi
    
    # Copy the entire rds folder from INPUT_FOLDER to SUBDIR
    if [ -d "$INPUT_FOLDER/rds" ]; then
        cp -r "$INPUT_FOLDER/rds" "$SUBDIR"
        
        # Inside the copied directory, delete any file containing "MAUDE"
        find "$SUBDIR/rds" -type f -name "*MAUDE*" -exec rm {} \;
    fi
    
    # Create QC_filtered_subsampled directory
    QC_DIR="$SUBDIR/QC_filtered"
    mkdir -p "$QC_DIR"
    now="$(date +"%T")"
    echo "$now   Subdirectories created successfully!"
    
    # Subsample .fastq.gz or .txt.gz files in QC_filtered directory
    if [ -d "$INPUT_FOLDER/QC_filtered" ]; then
        for file in "$INPUT_FOLDER/QC_filtered"/*; do
            # Check for .fastq.gz files and subsample them
            if [[ "$file" == *.fastq.gz ]]; then

                # Use seqtk to randomly subsample the file and save it in the subsampled directory
                seqtk sample "$file" "$percent_decimal" | gzip > "$QC_DIR/$(basename "$file")"
                now="$(date +"%T")"
                echo "$now   Percent: $percent, Decimal: $percent_decimal"
            # Check for .txt.gz files and subsample them
                now="$(date +"%T")"
                echo "$now   Subsampeling of file $file completed"
            else
                now="$(date +"%T")"
                echo "$now   $file is not a .fastq.gz file and will not be subsampled"
            fi
        done
    fi
    now="$(date +"%T")"
    echo "$now   subsampled files created successfully!"
done
now="$(date +"%T")"
echo "$now   All Done!"

