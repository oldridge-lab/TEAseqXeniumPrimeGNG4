#!/bin/bash

#SBATCH --job-name=xp_project_sketch             # Job name
#SBATCH --output=log/xp_project_sketch_log.txt   # Log file
#SBATCH --error=log/xp_project_sketch_status.txt # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                      # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=16                       # Number of CPU cores per task
#SBATCH --mem=512G                               # Memory per node
#SBATCH --time=0-12:00:00                        # Wall clock time (D-HH:MM:SS)

cd /filepath/xenium_data_processing/step6_project_sketch/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Running script"
Rscript step6_project_sketch.R

echo "Bash script complete"