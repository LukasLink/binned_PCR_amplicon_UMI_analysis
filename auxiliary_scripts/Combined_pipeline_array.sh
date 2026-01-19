#!/bin/bash
#SBATCH -J CP_array_1     # Job Name                # chmod +x ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/Combined_pipeline_array.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/Combined_pipeline_array.sh
#SBATCH --mem 128g               # Total memory required for the job #  --dependency=afterok:46151850
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 12                   # Number of CPUs
#SBATCH -t 02:15:00             # Runtime until the job is forcefully canceled
#SBATCH --array=0-8
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/CP_array_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/CP_array_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

# Array Options
################################################################################
percentages=("90" "80" "70" "60" "50" "40" "30" "20" "10")

i=${SLURM_ARRAY_TASK_ID}          # typically 0..8 if you have 9 entries
pct=${percentages[$i]}            # e.g. "90"

OUTPUT_FOLDER_ROOT="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/1/HepG2_dual_rep_PA_subsample"
OUTPUT_FOLDER="${OUTPUT_FOLDER_ROOT}/subsample_${pct}"

echo "Task $i -> ${OUTPUT_FOLDER}"
################################################################################

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################

# Input and output folder paths
INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/RAW/Liangfu_iBeer_2/Atto/"
# OUTPUT_FOLDER_ROOT="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/3/HepG2_dual_rep_PA_subsample/subsample_10"

# Regex pattern for extracting UMIs
#REGEX_PATTERN="^(?P<discard_1>.{0,5})(?P<umi_1>.{10})(?P<discard_2>[AGC]{2})(?P<discard_3>GTGGAAAGGACGAAACACCG){e<=1}"
# Regex pattern for extracting Reads
REGEX_PATTERN="^(?P<discard_1>.{0,4})(?P<umi_1>.{1})(?P<discard_2>TCTTGTGGAAAGGACGAAACACCG){e<=1}"

# UMI-tools needs the UMIs to be removed from the reads and added to the Readname,s
# it does this by capturing the UMIs as "umi_1", "umi_2"" etc. and throws away 
# everything captured by "discard_1", "discard_2", etc. 

# STAR settings
SJDB_OVERHANG=56 #This is the lenght read -1
Genome_SA_index_N_Bases=10 #Calculated by Combined_pipline_support.rmd
NUM_THREADS=10 # The number of Threads used by star, should NOT be higher than n set above. 
MAX_MEM=100000000000 # The number of bytes available to STAR, should NOT be higher than mem set above

# Decide if we process UMIs (leave false) or reads (set to true)
READS=true #if processing Reads, set all except Grouping, Filtering, deduplication to false

# Skip Options
SKIP_QC=true
SKIP_UMI_EXTRACTION=true
SKIP_GENOME_GENERATE=false
SKIP_MAPPING=false
SKIP_GROUPING=true
SKIP_FILTERING=true
SKIP_DEDUPLICATION=true
SKIP_IDXSTATS=false

# Decide if intermediary files should be removed to free up storage space
# particularly the large .tsv files can be removed if no read filtering is to be done. 
CLEAN_INTERMEDIARIES=false
CLEAN_TSVs=false # will automatically be set to true if CLEAN_INTERMEDIARIES is true. 

################################################################################
# END OF USER OPTIONS
################################################################################

# Function to adjust the paths to make them "//" proof
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
INPUT_FOLDER="$(normalize_path "$INPUT_FOLDER")"
OUTPUT_FOLDER="$(normalize_path "$OUTPUT_FOLDER")"

# Define necessary subfolders
QC_FILTERED="$OUTPUT_FOLDER/QC_filtered"
UMI_EXTRACTED="$OUTPUT_FOLDER/UMI_extracted"
STAR_INDEX="$OUTPUT_FOLDER/star_index"
MAPPED="$OUTPUT_FOLDER/mapped"
GROUPED="$OUTPUT_FOLDER/grouped"
READ_FILTERED="$OUTPUT_FOLDER/Read_filtered"
DEDUP="$OUTPUT_FOLDER/dedup"
LOGS="$OUTPUT_FOLDER/logs"

# Create necessary output directories
mkdir -p "$QC_FILTERED" "$UMI_EXTRACTED" "$STAR_INDEX" "$MAPPED" "$GROUPED" "$READ_FILTERED" "$DEDUP" "$LOGS" "$STAR_INDEX/NOPE"

# Rename .txt.gz files to .fastq.gz in INPUT_FOLDER
for f in "$INPUT_FOLDER"/*.txt.gz; do
    [[ -e "$f" ]] || continue  # skip if no matches
    new_name="${f%.txt.gz}.fastq.gz"
    now="$(date +"%T")"
    echo "$now   Renaming: $f -> $new_name"
    mv -- "$f" "$new_name"
done

if [ "$SKIP_QC" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting QC"
    # Load modules
    module purge
    ml seqtk/1.3-GCC-11.2.0
    
    # Process each file in the input folder
    for file in "$INPUT_FOLDER"/*.fastq.gz; do
        filename=$(basename "$file")
        sample_name="${filename%.fastq.gz}"
        seqtk seq -q20 -Q20 -L"$SJDB_OVERHANG" -n N "$file" | gzip > "$QC_FILTERED/$filename"
    done
    now="$(date +"%T")"
    echo "$now   Finished QC"
else
    now="$(date +"%T")"
    echo "$now   Skipping QC step as SKIP_QC is set to true"
fi

if [ "$SKIP_UMI_EXTRACTION" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting UMI Extraction"
    # Load UMI-tools module
    module purge
    ml UMI-tools/1.1.2-foss-2021b-Python-3.9.6
    
    
    # Extract UMIs
    for file in "$QC_FILTERED"/*.fastq.gz; do
        filename=$(basename "$file")
        sample_name="${filename%.fastq.gz}"
        umi_tools extract --stdin "$file" --stdout "$UMI_EXTRACTED/$filename" --extract-method=regex --bc-pattern="$REGEX_PATTERN"
    done
    now="$(date +"%T")"
    echo "$now   Finished UMI extraction"
else
    now="$(date +"%T")"
    echo "$now   Skipping UMI_EXTRACTION step as SKIP_UMI_EXTRACTION is set to true"
fi

if [ "$READS" == "true" ]; then
  MAPPING_SOURCE_DIR="$QC_FILTERED"
  now="$(date +"%T")"
  echo "$now   Using QC_filtered for mapping"
else
  MAPPING_SOURCE_DIR="$UMI_EXTRACTED"
  now="$(date +"%T")"
  echo "$now   Using UMI_extracted for mapping"
fi

if [ "$SKIP_GENOME_GENERATE" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting Genome Generation"
    # Load STAR module
    module purge
    ml STAR/2.7.11b-GCC-13.2.0
    
    # Generate genome indices
    STAR --runThreadN "$NUM_THREADS" --runMode genomeGenerate --genomeDir "$STAR_INDEX/NOPE" \
         --genomeFastaFiles "$OUTPUT_FOLDER/ref/NOPE_ref.fa" \
         --sjdbGTFfile "$OUTPUT_FOLDER/ref/NOPE_ref.gtf" \
         --limitGenomeGenerateRAM "$MAX_MEM" \
         --genomeSAindexNbases "$Genome_SA_index_N_Bases" --sjdbOverhang "$SJDB_OVERHANG"
    # IMPORTANT:
    # For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical
    # value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal
    # to 9, for 100 kiloBase genome, this is equal to 7.
    
    # and --sjdbOverhang should be ReadLength-1
    now="$(date +"%T")"
    echo "$now   Finished Genome Generation"
else
    now="$(date +"%T")"
    echo "$now   Skipping GENOME_GENERATE step as SKIP_GENOME_GENERATE is set to true"
fi

if [ "$SKIP_MAPPING" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting Mapping"
    # Map reads
    # Load STAR module
    module purge
    ml STAR/2.7.11b-GCC-13.2.0
    for file in "$MAPPING_SOURCE_DIR"/*.fastq.gz; do
        filename=$(basename "$file")
        sample_name="${filename%.fastq.gz}"

        # Build STAR args depending on SJDB_OVERHANG
        star_extra_args=()
    
        if (( SJDB_OVERHANG >= 70 )); then
            # Case 1: use your original snippet (no extra filters)
            :
        elif (( SJDB_OVERHANG >= 50 )); then
            # Case 2: 70 > SJDB_OVERHANG >= 50  (read_len 51–70)
            star_extra_args+=(
                --outFilterMatchNmin 30
                --outFilterScoreMinOverLread 0.5
                --outFilterMatchNminOverLread 0.5
            )
        elif (( SJDB_OVERHANG >= 30 )); then
            # Case 3: 50 > SJDB_OVERHANG >= 30  (read_len 31–50)
            star_extra_args+=(
                --outFilterMatchNmin 20
                --outFilterScoreMinOverLread 0.3
                --outFilterMatchNminOverLread 0.3
            )
        else
            # Case 4: SJDB_OVERHANG < 30  (read_len <= 30) -> hard stop
            now="$(date +"%T")"
            echo "$now   ERROR: Reads too short for STAR (SJDB_OVERHANG=${SJDB_OVERHANG}) in file: $file" >&2
            exit 1
        fi
    
        STAR --runThreadN "$NUM_THREADS" \
             --genomeDir "$STAR_INDEX/NOPE" \
             --readFilesCommand zcat \
             --readFilesIn "$file" \
             --outFileNamePrefix "$MAPPED/${sample_name}_" \
             --outSAMtype BAM SortedByCoordinate \
             "${star_extra_args[@]}"
    done
    now="$(date +"%T")"
    echo "$now   Finished Mapping"
    # Load SAMtools
    module purge
    ml SAMtools/1.21-GCC-13.3.0
    now="$(date +"%T")"
    echo "$now   Starting indexing the mapped bam files"
    # Index BAM files
    for file in "$MAPPED"/*.bam; do
        samtools index "$file"
    done
    now="$(date +"%T")"
    echo "$now   Finished indexing the mapped bam files"
else
    now="$(date +"%T")"
    echo "$now   Skipping MAPPING step as SKIP_MAPPING is set to true"
fi

if [ "$SKIP_GROUPING" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting Grouping (this will create some big .tsv files"
    # Group reads
    module purge
    ml UMI-tools/1.1.2-foss-2021b-Python-3.9.6
    
    for file in "$MAPPED"/*.bam; do
        filename=$(basename "$file")
        sample_name="${filename%_Aligned.sortedByCoord.out.bam}"
        umi_tools group -I "$file" --method adjacency --log "$GROUPED/$sample_name.log" \
                       --group-out "$GROUPED/$sample_name.tsv" --output-bam -S "$GROUPED/$sample_name.bam"
    done
    now="$(date +"%T")"
    echo "$now   Finished Grouping"
    # Index grouped BAM files
    module purge
    ml SAMtools/1.21-GCC-13.3.0
    for file in "$GROUPED"/*.bam; do
        samtools index "$file"
    done
else
    now="$(date +"%T")"
    echo "$now   Skipping GROUPING step as SKIP_GROUPING is set to true"
fi

### IMPORTANT!
# run Filter_UMIs_with_few_reads before proceeding, to generate _treshold_umis.txt

# If filtering is not skipped, proceed with filtering
if [ "$SKIP_FILTERING" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting Filtering based on custom threshold.tsv"
    # Load Picard module
    module purge
    ml picard/3.1.0-Java-17

    ### IMPORTANT!
    # run Filter_UMIs_with_few_reads before proceeding, to generate _treshold_umis.txt

    # Filter reads
    for file in "$GROUPED"/*.bam; do
        filename=$(basename "$file")
        sample_name="${filename%.bam}"
        java -jar $EBROOTPICARD/picard.jar FilterSamReads -I "$file" -O "$READ_FILTERED/$filename" \
               --READ_LIST_FILE "$GROUPED/${sample_name}_threshold_umis.txt" --FILTER excludeReadList
    done

    # Index filtered BAM files
    module purge
    ml SAMtools/1.21-GCC-13.3.0
    for file in "$READ_FILTERED"/*.bam; do
        samtools index "$file"
    done

    # Set source directory for deduplication
    DEDUP_SOURCE_DIR="$READ_FILTERED"
    now="$(date +"%T")"
    echo "$now   Finished Read_filtering"
else
    # If skipping filtering, use grouped BAMs directly
    DEDUP_SOURCE_DIR="$GROUPED"
    now="$(date +"%T")"
    echo "$now   Skipping Read_filtering step as SKIP_READ_FILTER is set to true"
fi

if [ "$SKIP_DEDUPLICATION" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting Deduplication"
    # Deduplicate reads
    module purge
    ml UMI-tools/1.1.2-foss-2021b-Python-3.9.6
    for file in "$DEDUP_SOURCE_DIR"/*.bam; do
        filename=$(basename "$file")
        sample_name="${filename%.bam}"
        umi_tools dedup --stdin="$file" --log="$DEDUP/$sample_name.log" --output-stats="$DEDUP/$sample_name" \
                       --method adjacency -S "$DEDUP/${sample_name}_dedup.bam"
    done
    now="$(date +"%T")"
    echo "$now   Finished Deduplication"
else
    now="$(date +"%T")"
    echo "$now   Skipping DEDUPLICATION step as SKIP_DEDUPLICATION is set to true"
fi

if [ "$SKIP_IDXSTATS" != "true" ]; then
    now="$(date +"%T")"
    echo "$now   Starting Indexing"

    if [ "$READS" == "true" ]; then
      INDEX_SOURCE_DIR="$MAPPED"
    else
      INDEX_SOURCE_DIR="$DEDUP"
    fi
    # Index deduplicated BAM files
    module purge
    ml SAMtools/1.21-GCC-13.3.0
    for file in "$INDEX_SOURCE_DIR"/*.bam; do
      filename=$(basename "$file" .bam)
      samtools idxstats "$file" > "$INDEX_SOURCE_DIR/${filename}_idxstats.txt"
    done
    now="$(date +"%T")"
    echo "$now   Finished IDXSTATS"
else
    now="$(date +"%T")"
    echo "$now   Skipping IDXSTATS step as SKIP_IDXSTATS is set to true"
fi

if [ "$CLEAN_INTERMEDIARIES" == "true" ]; then
    now="$(date +"%T")"
    echo "$now   Removing all intermediary files. This is not recommended unless you are certain everything will work first try!"
    CLEAN_TSVs=true
    rm -r "$QC_FILTERED" "$UMI_EXTRACTED" "$STAR_INDEX" "$GROUPED" "$READ_FILTERED"
    
    if [ "$READS" == "true" ]; then
        rm -r "$DEDUP"
    else
        rm -r "$MAPPED"
    fi
    now="$(date +"%T")"    
    echo "$now   Finished removing all intermediary files."
fi
if [ "$CLEAN_TSVs" == "true" ]; then
    now="$(date +"%T")"
    echo "$now   Removing the large .tsv files from the grouped folder. They are needed if you plan to do reads per UMI analysis!"

    for file in "$GROUPED"/*.tsv; do
        rm "$file"
    done
    now="$(date +"%T")"    
    echo "$now   Finished removing the large .tsv files."
fi
now="$(date +"%T")"
echo "$now   Processing complete."