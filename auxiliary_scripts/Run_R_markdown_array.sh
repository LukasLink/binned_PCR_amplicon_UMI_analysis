#!/bin/bash
#SBATCH -J Run_R_markdown_array     # Job Name                # chmod +x /home/link/Run_R_markdown_array.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/Run_R_markdown_array.sh   # sbatch --dependency=afterok:24467230 /home/link/Run_R_markdown.sh
#SBATCH --mem 6g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 00:30:00             # Runtime until the job is forcefully canceled
#SBATCH --array=0-11
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/Run_R_markdown_array_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/Run_R_markdown_array_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

# Output format: choose "script", "pdf", or "html"
MODE="script"  # <-- change this to "pdf" or "html" or "script"as needed
# Path to the Rmd file
RMD_FILE="/home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/data_analysis_with_MAUDE.Rmd"
# RMD_FILE="/home/link/NB_EXP030/Combined_pipline_support.Rmd"

# Define combinations
pipelines=("lukas" "lukas" "lukas")
data_type=("umis" "umis" "umis" "umis" "umis" "umis" "reads" "reads" "reads" "reads" "reads" "reads")
methods=("sum" "sum" "rep" "rep" "" "" "sum" "sum" "rep" "rep" "" "")
norm_method=("control_median" "" "control_median" "" "control_median" "" "control_median" "" "control_median" "" "control_median" "")


# Get combo based on SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
# PIPELINE="${pipelines[$i]}"
METHOD="${methods[$i]}"
DATA_TYPE="${data_type[$i]}"
# EXTRA_SUFFIX="${extra_suffix[$i]}"
NORM_METHOD="${norm_method[$i]}"

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
  Rscript "$R_SCRIPT" --pipeline "lukas" --data_type "$DATA_TYPE" --method "$METHOD" --norm_method "$NORM_METHOD" --drop_0s FALSE --recover_input TRUE # --extra_suffix "$EXTRA_SUFFIX"


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