#!/bin/bash
#SBATCH -J RRMA_3     # Job Name                # chmod +x /home/link/Run_R_markdown_array.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/Run_R_markdown_array.sh   # sbatch --dependency=afterok:24467230 /home/link/Run_R_markdown.sh
#SBATCH --mem 6g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 00:30:00             # Runtime until the job is forcefully canceled
#SBATCH --array=0-8
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/RRMA_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/RRMA_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

################################################################################
# USER OPTIONS - Adjust these as needed
################################################################################

# Output format: choose "script", "pdf", or "html"
MODE="script"  # <-- change this to "pdf" or "html" or "script"as needed
# Choose if options or folders are to be alternated.
DIFFERENT_OR_SAME="same" # can be "same" or "different"
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
  
fi

if [ "$DIFFERENT_OR_SAME" == "same" ]; then
  echo "Running the same combination for multiple folders"
  
  percentages=("90" "80" "70" "60" "50" "40" "30" "20" "10")
  pct=${percentages[$i]}            # e.g. "90"
  
  OUTPUT_FOLDER_ROOT="/g/steinmetz/link/Amplicon_barcode_analysis/PA_subsampeling/3/HepG2_dual_rep_PA_subsample"
  OUTPUT_FOLDER="${OUTPUT_FOLDER_ROOT}/subsample_${pct}"
  
  PIPELINE="lukas"
  METHOD=""
  DATA_TYPE="reads"
  NORM_METHOD="control_median"
  EXTRA_SUFFIX=""
  
  echo "Task $i -> ${OUTPUT_FOLDER}"
fi
################################################################################
# END OF USER OPTIONS
################################################################################

module purge
ml R-bundle-Bioconductor
ml R-bundle-CRAN
ml Pandoc
ml texlive
ml ICU

export R_LIBS_USER="/home/link/R/x86_64-pc-linux-gnu-library/4.4"
export R_LIBS_USER="/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.4.2/MAUDE"

# Make sure temp dir exists
mkdir -p /scratch/link/temp

# Where to extract the R script version
R_SCRIPT="/scratch/link/temp/$(basename "$RMD_FILE" .Rmd)_${i}.R"

# Logic to handle different modes
if [ "$MODE" == "script" ]; then

  echo "Running Rmd as plain R script..."
  Rscript -e "knitr::purl('$RMD_FILE', output='$R_SCRIPT', documentation = 0)"
  Rscript "$R_SCRIPT" --output_folder "$OUTPUT_FOLDER" --pipeline "$PIPELINE" --data_type "$DATA_TYPE" --method "$METHOD" --norm_method "$NORM_METHOD" --drop_0s FALSE --recover_input TRUE # --extra_suffix "$EXTRA_SUFFIX"


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