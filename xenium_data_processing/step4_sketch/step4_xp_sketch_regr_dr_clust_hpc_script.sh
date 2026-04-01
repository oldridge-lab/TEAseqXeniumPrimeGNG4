#!/bin/bash

#SBATCH --job-name=xp_sketch_regr_dr_clust             # Job name
#SBATCH --output=log/xp_sketch_regr_dr_clust_log.txt   # Log file
#SBATCH --error=log/xp_sketch_regr_dr_clust_status.txt # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                            # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=16                             # Number of CPU cores per task
#SBATCH --mem=512G                                     # Memory per node
#SBATCH --time=0-12:00:00                              # Wall clock time (D-HH:MM:SS)

cd /filepath/xenium_data_processing/step4_sketch/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Running script"
Rscript step4_xp_sketch_regr_dr_clust.R

echo "Bash script complete"