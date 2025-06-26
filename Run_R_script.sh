#!/bin/bash
#SBATCH -J Run_R_script_filter_count     # Job Name                # chmod +x /home/link/Run_R_script.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/Run_R_script.sh   # sbatch --dependency=afterok:23764729 /home/link/Run_R_script.sh
#SBATCH --mem 64g               # Total memory required for the job
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 1                   # Number of CPUs
#SBATCH -t 3-00:00:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/Run_R_script_filter_count.out
#SBATCH -e /g/steinmetz/link/logs/Run_R_script_filter_count.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail


# Path to the .R file
R_FILE="/home/link/NB_EXP030/filter_umis_with_few_reads.R"

module purge
ml R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1

Rscript "$R_FILE" 

$1 \
$2 \
$3 \
$4