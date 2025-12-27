#!/bin/bash

#SBATCH --job-name=tc656a_qc                  # Job name
#SBATCH --output=log/tc656a_qc_log.txt        # Log file
#SBATCH --error=log/tc656a_qc_status.txt      # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                   # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=8                     # Number of CPU cores per task
#SBATCH --mem=128G                            # Memory per node
#SBATCH --time=0-04:00:00                     # Wall clock time (D-HH:MM:SS)

cd /filepath/xenium_data_processing/step2_qc/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Running script"
Rscript tc656a_qc.R

echo "Bash script complete"