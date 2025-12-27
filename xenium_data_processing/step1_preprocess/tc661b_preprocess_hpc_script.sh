#!/bin/bash

#SBATCH --job-name=tc661b_preprocess_script   # Job name
#SBATCH --output=log/tc661b_preproc_log.txt   # Log file
#SBATCH --error=log/tc661b_preproc_status.txt # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                   # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=16                    # Number of CPU cores per task
#SBATCH --mem=256G                            # Memory per node
#SBATCH --time=0-04:00:00                     # Wall clock time (D-HH:MM:SS)

cd /filepath/xenium_data_processing/step1_preprocess/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Running script"
Rscript tc661b_preprocess.R

echo "Bash script complete"