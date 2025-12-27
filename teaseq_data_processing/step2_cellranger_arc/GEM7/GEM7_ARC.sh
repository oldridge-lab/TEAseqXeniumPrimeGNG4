#!/bin/bash
#SBATCH --job-name=GEM7_ARC
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
#SBATCH --output=log/cr_log_%j.txt

echo "starting cellranger-arc count"

export export PATH=/mnt/isilon/oldridge_lab/install/cellranger-arc-2.0.2:$PATH

cd /output_folder_filepath/

cellranger-arc count --id=GEM7_ARC\
        --reference=/reference_genome_fliepath/refdata-cellranger-arc-GRCh38-2020-A-2.0.0\
        --libraries=/flowcell_output_filepath/GEM7_ARC_libraries.csv\
        
echo "done"