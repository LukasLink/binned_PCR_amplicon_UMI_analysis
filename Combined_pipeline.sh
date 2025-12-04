#!/bin/bash
#SBATCH -J Combined_pipeline_split     # Job Name                # chmod +x /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/Combined_pipeline.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/Combined_pipeline.sh
#SBATCH --mem 128g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 12                   # Number of CPUs
#SBATCH -t 05:00:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/log_Combined_pipeline_split.out
#SBATCH -e /g/steinmetz/link/logs//log_Combined_pipeline_split.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################

# Input and output folder paths
INPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/fake_untreated"
OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/fake_untreated/"
# Regex pattern for extracting UMIs
# REGEX_PATTERN="^(?P<discard_1>.{0,5})(?P<umi_1>.{10})(?P<discard_2>[AGC]{2})(?P<discard_3>GTGGAAAGGACGAAACACCG){e<=1}"
# Regex pattern for extracting Reads
REGEX_PATTERN="^(?P<discard_1>.{0,4})(?P<umi_1>.{1})(?P<discard_2>TCTTGTGGAAAGGACGAAACACCG){e<=1}"
# UMI-tools needs the UMIs to be removed from the reads and added to the Readname,
# it does this by capturing the UMIs as "umi_1", "umi_2"" etc. and throws away 
# everything captured by "discard_1", "discard_2", etc. 

# STAR settings
SJDB_OVERHANG=121 #This is the lenght read -1
Genome_SA_index_N_Bases=9 #Calculated by Combined_pipline_support.rmd
NUM_THREADS=10 # The number of Threads used by star, should NOT be higher than n set above. 
MAX_MEM=100000000000 # The number of bytes available to STAR, should NOT be higher than mem set above

# Skip Options
# Decide if we want to skip read filtering or not (if yes will exclude UMIs below a certain number of reads)
READS=true #if processing Reads, set all except Grouping, Filtering, deduplication to false
SKIP_QC=true
SKIP_UMI_EXTRACTION=false
SKIP_GENOME_GENERATE=false
SKIP_MAPPING=false
SKIP_GROUPING=true
SKIP_FILTERING=true
SKIP_DEDUPLICATION=true
SKIP_IDXSTATS=false
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
    echo "Renaming: $f -> $new_name"
    mv -- "$f" "$new_name"
done

if [ "$SKIP_QC" != "true" ]; then
    # Load modules
    module purge
    ml seqtk/1.3-GCC-11.2.0
    
    # Process each file in the input folder
    for file in "$INPUT_FOLDER"/*.fastq.gz; do
        filename=$(basename "$file")
        sample_name="${filename%.fastq.gz}"
        seqtk seq -q20 -Q20 -L100 -n N "$file" | gzip > "$QC_FILTERED/$filename"
    done
    echo "Finished QC"
else
    echo "Skipping QC step as SKIP_QC is set to true"
fi

if [ "$SKIP_UMI_EXTRACTION" != "true" ]; then
    # Load UMI-tools module
    module purge
    ml UMI-tools/1.1.2-foss-2021b-Python-3.9.6
    
    
    # Extract UMIs
    for file in "$QC_FILTERED"/*.fastq.gz; do
        filename=$(basename "$file")
        sample_name="${filename%.fastq.gz}"
        umi_tools extract --stdin "$file" --stdout "$UMI_EXTRACTED/$filename" --extract-method=regex --bc-pattern="$REGEX_PATTERN"
    done
    echo "Finished UMI extraction"
else
    echo "Skipping UMI_EXTRACTION step as SKIP_UMI_EXTRACTION is set to true"
fi

if [ "$SKIP_GENOME_GENERATE" != "true" ]; then
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
    echo "Finished Genome Generation"
else
    echo "Skipping GENOME_GENERATE step as SKIP_GENOME_GENERATE is set to true"
fi

if [ "$SKIP_MAPPING" != "true" ]; then
    # Map reads
    # Load STAR module
    module purge
    ml STAR/2.7.11b-GCC-13.2.0
    for file in "$UMI_EXTRACTED"/*.fastq.gz; do
        filename=$(basename "$file")
        sample_name="${filename%.fastq.gz}"
        STAR --runThreadN "$NUM_THREADS" --genomeDir "$STAR_INDEX/NOPE" --readFilesCommand zcat \
             --readFilesIn "$file" --outFileNamePrefix "$MAPPED/${sample_name}_" --outSAMtype BAM SortedByCoordinate
    done
    echo "Finished Mapping"
    # Load SAMtools
    module purge
    ml SAMtools/1.21-GCC-13.3.0
    
    # Index BAM files
    for file in "$MAPPED"/*.bam; do
        samtools index "$file"
    done
else
    echo "Skipping MAPPING step as SKIP_MAPPING is set to true"
fi

if [ "$SKIP_GROUPING" != "true" ]; then
  # Group reads
  module purge
  ml UMI-tools/1.1.2-foss-2021b-Python-3.9.6
  
  for file in "$MAPPED"/*.bam; do
      filename=$(basename "$file")
      sample_name="${filename%_Aligned.sortedByCoord.out.bam}"
      umi_tools group -I "$file" --method adjacency --log "$GROUPED/$sample_name.log" \
                     --group-out "$GROUPED/$sample_name.tsv" --output-bam -S "$GROUPED/$sample_name.bam"
  done
  echo "Finished Grouping"
  # Index grouped BAM files
  module purge
  ml SAMtools/1.21-GCC-13.3.0
  for file in "$GROUPED"/*.bam; do
      samtools index "$file"
  done
else
    echo "Skipping GROUPING step as SKIP_GROUPING is set to true"
fi

### IMPORTANT!
# run Filter_UMIs_with_few_reads before proceeding, to generate _treshold_umis.txt

# If filtering is not skipped, proceed with filtering
if [ "$SKIP_FILTERING" != "true" ]; then
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
    echo "Finished Read_filtering"
else
    # If skipping filtering, use grouped BAMs directly
    DEDUP_SOURCE_DIR="$GROUPED"
    echo "Skipping Read_filtering step as SKIP_READ_FILTER is set to true"
fi

if [ "$SKIP_DEDUPLICATION" != "true" ]; then
  # Deduplicate reads
  module purge
  ml UMI-tools/1.1.2-foss-2021b-Python-3.9.6
  for file in "$DEDUP_SOURCE_DIR"/*.bam; do
      filename=$(basename "$file")
      sample_name="${filename%.bam}"
      umi_tools dedup --stdin="$file" --log="$DEDUP/$sample_name.log" --output-stats="$DEDUP/$sample_name" \
                     --method adjacency -S "$DEDUP/${sample_name}_dedup.bam"
  done
  echo "Finished Deduplication"
else
    echo "Skipping DEDUPLICATION step as SKIP_DEDUPLICATION is set to true"
fi

if [ "$SKIP_IDXSTATS" != "true" ]; then

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
  echo "Finished IDXSTATS"
else
    echo "Skipping IDXSTATS step as SKIP_IDXSTATS is set to true"
fi
echo "Processing complete."