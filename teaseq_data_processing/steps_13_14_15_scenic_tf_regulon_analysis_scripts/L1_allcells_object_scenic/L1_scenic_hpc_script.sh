#!/bin/bash
#SBATCH --job-name=L1_teaseq_obj_scenic              # Job name
#SBATCH --output=log/L1_teaseq_obj_scenic_log.txt    # Log file
#SBATCH --error=log/L1_teaseq_obj_scenic_status.txt  # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                          # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=800G                                   # Memory per node
#SBATCH --time=4-00:00:00                            # Wall clock time (D-HH:MM:SS)

# R SCENIC Script for L1 TEAseq object
echo "Starting script"

cd /filepath/step13_bulk_harmony/SCENIC/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Execute SCENIC R script"
Rscript L1_teaseq_obj_scenic_script.R

echo "Done"