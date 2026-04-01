#!/bin/bash
#SBATCH --job-name=L3_teaseq_obj_scenic              # Job name
#SBATCH --output=log/L3_teaseq_obj_scenic_log.txt    # Log file
#SBATCH --error=log/L3_teaseq_obj_scenic_status.txt  # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                          # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=20                           # Number of CPU cores per task
#SBATCH --mem=256G                                   # Memory per node
#SBATCH --time=0-48:00:00                            # Wall clock time (D-HH:MM:SS)

# R SCENIC Script for L3 TEAseq object
echo "Starting script"

cd /filepath/step15_tfh_subcluster/SCENIC/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Execute SCENIC R script"
Rscript L3_teaseq_obj_scenic_script.R

echo "Done"