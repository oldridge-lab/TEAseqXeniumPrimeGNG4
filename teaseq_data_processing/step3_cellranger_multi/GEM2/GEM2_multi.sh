#!/bin/bash
#SBATCH --job-name=GEM2_multi
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --time=8:00:00
#SBATCH --output=log/cr_log_%j.txt

echo "starting cellranger multi"

export PATH=/install_path/cellranger-8.0.0:$PATH

cd /output_filepath/

cellranger multi --id=GEM2_multi\
        --csv=/config_spreadsheet_filepath/GEM2_multi_config.csv\
        
echo "done"