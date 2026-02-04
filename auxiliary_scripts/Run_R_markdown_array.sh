#!/bin/bash
#SBATCH -J RRMA_rep_GALNAC     # Job Name                # chmod +x /home/link/Run_R_markdown_array.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/Run_R_markdown_array.sh   
#SBATCH --mem 6g               # Total memory required for the job    # sbatch --dependency=afterok:24467230 /home/link/Run_R_markdown.sh
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 00:08:00             # Runtime until the job is forcefully canceled
#SBATCH --array=0-18
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/%x_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################

# WARNING: I am currently hard coding th R-libs Path below the user options!!!

# Output format: choose "script", "pdf", or "html"
MODE="script"  # <-- change this to "pdf" or "html" or "script"as needed
# Choose if options or folders are to be alternated.
DIFFERENT_OR_SAME="same" # can be "same" or "different"
SAME_OPTIONS="replicates" # can be "replicates", "directories", or "subsample", "sublib_skip"
# Path to the Rmd file
RMD_FILE="/home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/data_analysis_with_MAUDE.Rmd"

i=$SLURM_ARRAY_TASK_ID
# For iterating through different parameter combos
if [ "$DIFFERENT_OR_SAME" == "different" ]; then
  echo "Running multiple different combinations"
  
  # Define combinations must all be the same length
  pipelines=("lukas" "lukas" "lukas")
  data_type=("umis" "umis" "umis" "umis" "umis" "umis" "reads" "reads" "reads" "reads" "reads" "reads")
  methods=("sum" "sum" "rep" "rep" "" "" "sum" "sum" "rep" "rep" "" "")
  norm_method=("control_median" "" "control_median" "" "control_median" "" "control_median" "" "control_median" "" "control_median" "")
  
  # Assign combination for this run
  OUTPUT_FOLDER="path/to/your/folder"
  PIPELINE="${pipelines[$i]}"
  METHOD="${methods[$i]}"
  DATA_TYPE="${data_type[$i]}"
  NORM_METHOD="${norm_method[$i]}"
  EXTRA_SUFFIX=""
  SKIP_SUBLIB=""
  
fi

if [ "$DIFFERENT_OR_SAME" == "same" ]; then
  echo "Running the same combination for multiple folders"
  
  if [ "$SAME_OPTIONS" == "subsample" ]; then
    percentages=("90" "80" "70" "60" "50" "40" "30" "20" "10")
    pct=${percentages[$i]}            # e.g. "90"
    
    OUTPUT_FOLDER_ROOT="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/3/HepG2_dual_rep_PA_subsample"
    OUTPUT_FOLDER="${OUTPUT_FOLDER_ROOT}/subsample_${pct}"
    
    PIPELINE="lukas"
    METHOD=""
    DATA_TYPE="reads"
    NORM_METHOD="control_median"
    EXTRA_SUFFIX=""
    SKIP_SUBLIB=""
    
    echo "Susbsample Task $i -> ${OUTPUT_FOLDER}"
  fi
  
  if [ "$SAME_OPTIONS" == "replicates" ]; then
  
    # OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_GALNAC"
    OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_GALNAC"
    
    PIPELINE="lukas"
    METHOD=""
    DATA_TYPE="reads"
    NORM_METHOD="control_median"
    EXTRA_SUFFIX="rep$((i+1))"
    SKIP_SUBLIB=""
    
    echo "MAUDE Replicate Task $i -> $EXTRA_SUFFIX"
  fi
  
  if [ "$SAME_OPTIONS" == "directories" ]; then
  
    # OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_GALNAC"
    OUTPUT_FOLDERS=("/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_GALNAC" "/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_PA" "/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_DCA")
    OUTPUT_FOLDER=${OUTPUT_FOLDERS[$i]} 
    
    PIPELINE="lukas"
    METHOD=""
    DATA_TYPE="reads"
    NORM_METHOD="control_median"
    EXTRA_SUFFIX=""
    SKIP_SUBLIB=""
    
    echo "Directory Task $i -> $OUTPUT_FOLDER"
  fi
  
  if [ "$SAME_OPTIONS" == "sublib_skip" ]; then
  
    # OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_GALNAC"
    OUTPUT_FOLDER="/g/steinmetz/link/Amplicon_barcode_analysis/HepG2_dual_rep_GALNAC"
    SKIP_SUBLIB_LIST=("L1" "L2" "L3" "L4")
    SKIP_SUBLIB=${SKIP_SUBLIB_LIST[$i]}
    
    
    PIPELINE="lukas"
    METHOD=""
    DATA_TYPE="reads"
    NORM_METHOD="control_median"
    EXTRA_SUFFIX=""
    
    echo "Susbsample Task $i -> $SKIP_SUBLIB_LIST"
  fi
fi
################################################################################
# END OF USER OPTIONS
################################################################################

module purge
ml R/4.4.2-gfbf-2024a
# ml R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
# ml R-bundle-CRAN/2024.11-foss-2024a
# ml Pandoc
# ml texlive
# ml ICU

# Make sure the library paths are correct
unset R_LIBS
unset R_LIBS_SITE
# export R_LIBS_USER="/home/link/R/x86_64-pc-linux-gnu-library/4.4"
# export R_LIBS_USER="/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.4.2/MAUDE"
export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/glibc-2.28/R-4.4.2
export R_LIBS=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/glibc-2.28/R-4.4.2
export R_LIBS_SITE=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/glibc-2.28/R-4.4.2

echo "R_LIBS"
echo "$R_LIBS"
echo "R_LIBS_SITE"
echo "$R_LIBS_SITE"
echo "R_LIBS_USER"
echo "$R_LIBS_USER"
# Make sure temp dir exists
mkdir -p /scratch/link/temp

# Where to extract the R script version
R_SCRIPT="/scratch/link/temp/$(basename "$RMD_FILE" .Rmd)_${i}.R"

# Logic to handle different modes
if [ "$MODE" == "script" ]; then

  echo "Running Rmd as plain R script..."
  Rscript --vanilla -e "knitr::purl('$RMD_FILE', output='$R_SCRIPT', documentation = 0)"
  export SOURCE_RMD="$RMD_FILE"
  Rscript --vanilla "$R_SCRIPT" --first_time T --output_folder "$OUTPUT_FOLDER" --pipeline "$PIPELINE" --data_type "$DATA_TYPE" --method "$METHOD" --norm_method "$NORM_METHOD" --drop_0s FALSE --recover_input TRUE --extra_suffix "$EXTRA_SUFFIX" --skip_list_sublib "$SKIP_SUBLIB"


elif [ "$MODE" == "pdf" ]; then
  echo "Rendering Rmd to PDF..."
  Rscript -e "rmarkdown::render('$RMD_FILE', output_format = 'pdf_document')"

elif [ "$MODE" == "html" ]; then
  echo "Rendering Rmd to HTML..."
  Rscript -e "rmarkdown::render('$RMD_FILE', output_format = 'html_document')"

else
  echo "Invalid MODE: $MODE. Use 'script', 'pdf', or 'html'."
  exit 1
fi



$1 \
$2 \
$3 \
$4