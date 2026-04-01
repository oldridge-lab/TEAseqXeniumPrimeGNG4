#!/bin/bash

#SBATCH --job-name=snp_demux_gem6       # Job name
#SBATCH --output=snp_demux_gem6         # Output file
#SBATCH --error=snp_demux_gem6_error    # Error file
#SBATCH --ntasks-per-node=1             # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=32              # Number of CPU cores per task
#SBATCH --mem=50G                       # Memory per node
#SBATCH --time=1-00:00:00               # Wall clock time (D-HH:MM:SS)

cd /filepath/step4_snp_demultiplex/gem6

export PATH=/install_filepath/souporcell_latest.sif:$PATH

# Create unique output directory for selected GEM well
SOUPORCELL_OUTDIR="/filepath/step4_snp_demultiplex/gem6"
mkdir -p "$SOUPORCELL_OUTDIR"

# Provide path to corresponding unzipped cellranger-arc filtered cell barcode list for selected GEM well 
BARCODES="/filepath/step4_snp_demultiplex/gem_barcodes/gem6/barcodes.tsv"
          
# Provide path to binary alignment map file for selected GEM well 
BAM="/filepath/cellranger_arc_output/GEM6_ARC/outs/gex_possorted_bam.bam"

# Execute the souporcell_pipeline.py command for each GEM well. Provide -f path to reference GEX genome from cellranger-arc, as well as common SNP file
singularity exec --bind /filepath/ \
            /install_filepath/souporcell_latest.sif souporcell_pipeline.py \
            -i "$BAM" \
            -b "$BARCODES" \
            -f /reference_genome_filepath/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
            -t 32 \
            -o "$SOUPORCELL_OUTDIR" \
            -k 8 \
            --common_variants /1kgenomes_variant_reference_filepath/filtered_2p_1kgenomes_unchr.vcf \
            
echo "done"