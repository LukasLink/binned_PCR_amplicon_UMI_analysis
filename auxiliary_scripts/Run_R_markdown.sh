#!/bin/bash
#SBATCH -J Run_R_markdown_maude     # Job Name                # chmod +x /home/link/Run_R_markdown.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/Run_R_markdown.sh   # sbatch --dependency=afterok:24956406 /home/link/Run_R_markdown.sh
#SBATCH --mem 8g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 00:20:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/Run_R_markdown_maude.out
#SBATCH -e /g/steinmetz/link/logs/Run_R_markdown_maude.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

# Output format: choose "script", "pdf", or "html"
MODE="script"  # <-- change this to "pdf" or "html" or "script"as needed
# Path to the Rmd file
RMD_FILE="/home/link/NB_EXP030/Hit_heatmaps.Rmd"
# RMD_FILE="/home/link/NB_EXP030/NB_EXP030.Rmd"
# RMD_FILE="/home/link/NB_EXP030/Combined_pipline_support.Rmd"

module purge
ml R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1
ml R-bundle-CRAN/2024.06-foss-2023b
ml Pandoc
ml texlive/20240312-GCC-13.3.0
ml ICU/74.1-GCCcore-13.2.0

export R_LIBS_USER="/home/link/R/x86_64-pc-linux-gnu-library/4.4"
# Make sure temp dir exists
mkdir -p /scratch/link/temp

# Where to extract the R script version
R_SCRIPT="/scratch/link/temp/$(basename "$RMD_FILE" .Rmd).R"

# Logic to handle different modes
if [ "$MODE" == "script" ]; then
  echo "Running Rmd as plain R script..."
  Rscript -e "knitr::purl('$RMD_FILE', output='$R_SCRIPT', documentation = 0)"
  Rscript "$R_SCRIPT"
  # Rscript "$R_SCRIPT" --pipeline "lukas" --data_type "umis" --method "rep" --drop_0s FALSE --recover_input TRUE --extra_suffix ""

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