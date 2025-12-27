#!/bin/bash

#SBATCH --job-name=xp_donor_merge_script    # Job name
#SBATCH --output=log/merge_script_log.txt   # Log file
#SBATCH --error=log/merge_script_status.txt # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                 # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=16                  # Number of CPU cores per task
#SBATCH --mem=512G                          # Memory per node
#SBATCH --time=0-6:00:00                    # Wall clock time (D-HH:MM:SS)

cd /filepath/xenium_data_processing/step3_merge/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Running script"
Rscript step3_merge_xenium_donors.R

echo "Bash script complete"