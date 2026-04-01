#!/bin/bash

#SBATCH --job-name=step8b_xp_nncd4_subclust_project_sketch             # Job name
#SBATCH --output=log/step8b_xp_nncd4_subclust_project_sketch_log.txt   # Log file
#SBATCH --error=log/step8b_xp_nncd4_subclust_project_sketch_status.txt # Function status, warnings, and errors 
#SBATCH --ntasks-per-node=1                                            # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=16                                             # Number of CPU cores per task
#SBATCH --mem=400G                                                     # Memory per node
#SBATCH --time=0-6:00:00                                               # Wall clock time (D-HH:MM:SS)

cd /filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/

echo "Load R Version 4.4.0"
module load R/4.4.0

echo "Running R script"
Rscript step8b_xp_nncd4_subclust_project_sketch.R

echo "Bash script complete"