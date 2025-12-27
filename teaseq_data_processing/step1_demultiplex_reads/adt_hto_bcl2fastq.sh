#!/bin/bash
#SBATCH --job-name=adt_hto_bcl2fastq2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --output=/log/demux_log_%j.txt

echo "starting bcl2fastq"

cd /output_folder/

module load bcl2fastq2

bcl2fastq --use-bases-mask=Y28n*,I8n*,Y15n* \
--create-fastq-for-index-reads \
--minimum-trimmed-read-length=8 \
--mask-short-adapter-reads=8 \
--ignore-missing-positions \
--ignore-missing-filter \
--ignore-missing-bcls \
--barcode-mismatches=1 \
-r 24 -w 24 -p 80 \
-R /raw_reads_folder \
--output-dir /filepath/adt_hto_bcl2fastq_output \
--interop-dir /filepath/adt_hto_bcl2fastq_interop \
--sample-sheet /filepath/adt_hto_samplesheet.csv \
--no-lane-splitting \

echo "done"
