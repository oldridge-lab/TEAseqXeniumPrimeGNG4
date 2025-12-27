#!/bin/bash
#SBATCH --job-name=GEM1_ARC
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --output=log/cr_log_%j.txt

echo "starting cellranger-arc"

export export PATH=/mnt/isilon/oldridge_lab/install/cellranger-arc-2.0.2:$PATH

cd /output_folder_filepath/

cellranger-arc count --id=GEM1_ARC\
        --reference=/reference_genome_fliepath/refdata-cellranger-arc-GRCh38-2020-A-2.0.0\
        --libraries=/flowcell_output_filepath/GEM1_ARC_libraries.csv\
        
echo "done"