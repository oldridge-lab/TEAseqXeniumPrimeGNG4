# SCENIC analysis for Level 3 TEAseq Tfh-like cell object
# Code adapted from SCENIC tutorials (https://scenic.aertslab.org/tutorials/)

# 1) Set up R working environment ----

print ('Set up R working environment')

setwd('/filepath/step15_tfh_subcluster/SCENIC')

set.seed(26)

# Load packages and record version numbers
library(Matrix) # 1.7-0
library(tidyselect) # 1.2.1
library(BiocManager) # 1.30.23
library(GenomeInfoDb) # 1.40.0
library(GenomicRanges) # 1.56.0
library(IRanges) # 2.38.0
library(Rsamtools) # 2.20.0
library(S4Vectors) # 0.42.0
library(BiocGenerics) # 0.50.0
library(biovizBase) # 1.52.0
library(GenomicFeatures) # 1.56.0
library(EnsDb.Hsapiens.v86) # 2.99.0
library(BSgenome.Hsapiens.UCSC.hg38) # 1.4.5
library(Signac) # 1.13.0
library(Seurat) # 5.1.0
library(SeuratObject) # 5.0.2
library(data.table) # 1.15.4
library(readr) # 2.1.5
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4
library(stringr) # 1.5.1
library(glmGamPoi) # 1.16.0
library(TFBSTools) # 1.42.0
library(JASPAR2020) # 0.99.1
library(motifmatchr) # 1.26.0
library(ggseqlogo) # 0.2
library(caret) # 6.0-94
library(presto) # 1.0.0
library(clustree) # 0.5.1
library(SoupX) # 1.6.2
library(scDblFinder) # 1.19.1
library(Nebulosa) # 1.14.0
library(clusterProfiler) # 4.12.6 
library(enrichplot) # 1.24.2
library(org.Hs.eg.db) # 3.19.1
library(harmony) # 1.2.0
library(SCpubr) # 2.0.2
library(paletteer) # 1.6.0
library(RColorBrewer) # 1.1-3
library(EnhancedVolcano) # 1.22.0
library(AUCell) # 1.26.0
library(msigdbr) # 7.5.1
library(RcisTarget) # 1.23.1
library(GENIE3) # 1.26.0
library(SCENIC) # 1.3.1
library(arrow) # 18.1.0.1
library(doMC) # 1.3.8
library(R2HTML) # 2.3.4
library(reshape2) # 1.4.4

# 2) Get ranked gene x motif matrix from cisTarget database for SCENIC R implementation ----

print('Get ranked gene x motif matrix from cisTarget database for SCENIC R implementation')

# Download genes_vs_motifs.rankings.feather files from Aerts Lab hg38 refseq_r80 motif collection updated for SCENIC+
# URL (https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/)
# Matrix containing motifs as rows and genes as columns and ranking position for each gene and motif (based on CRM scores) as values. To be used with cisTarget (R).
# The search space around the TSS of the gene in which the motif is scored is indicated in the database name:
# TSS+/-10kb: 10kb around the TSS (total: 20kb).
# 500bpUp100Dw: 500bp upstream of TSS, and 100bp downstream.

# Need to adjust motifs column format of new cisTarget database files to format expected by SCENIC R implementation, then save modified feather files (from https://github.com/aertslab/SCENIC/issues/471)

# Load SCENIC feather file for TSS 10 kb up/downstream gene vs motif ranking - as downloaded from cisTarget database
tss.10kb.feather <- read_feather("/filepath/step15_tfh_subcluster/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
# Rename 'motifs' column to 'features' in dataframe
colnames(tss.10kb.feather)[colnames(tss.10kb.feather) == "motifs"] <- "features"
# Move 'features' column first in dataframe 
tss.10kb.feather <- tss.10kb.feather[, c("features", setdiff(colnames(tss.10kb.feather), "features"))]
# Save modified feather file for TSS 10 kb up/downstream gene vs motif ranking - now amenable to SCENIC R implementation
write_feather(tss.10kb.feather, "/filepath/step15_tfh_subcluster/SCENIC/modified_hg38_10kb_up_down.feather")

# Load SCENIC feather file for TSS 500bp upstream / 100 bp downstream gene vs motif ranking - as downloaded from cisTarget database
tss.500bp.100bp.feather <- read_feather("/filepath/step15_tfh_subcluster/SCENIC/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
# Rename 'motifs' column to 'features'
colnames(tss.500bp.100bp.feather)[colnames(tss.500bp.100bp.feather) == "motifs"] <- "features"
# Move 'features' column first in dataframe 
tss.500bp.100bp.feather <- tss.500bp.100bp.feather[, c("features", setdiff(colnames(tss.500bp.100bp.feather), "features"))]
# Save modified feather file for TSS 500bp upstream / 100 bp downstream gene vs motif ranking - now amenable to SCENIC R implementation
write_feather(tss.500bp.100bp.feather, "/filepath/step15_tfh_subcluster/SCENIC/modified_hg38_500bp_up_100bp_down.feather")

# 3) Import v10 motif annotations from RcisTarget ----
# From (https://www.bioconductor.org/packages/devel/bioc/manuals/RcisTarget/man/RcisTarget.pdf)
# motifAnnotations_hgnc: Annotations to HUMAN transcription factors for the rankings using motif collection version 10. Source: motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
data(list="motifAnnotations_hgnc", package="RcisTarget")

# Rename motif annotation object to avoid name mismatch error (from https://github.com/aertslab/SCENIC/issues/364)
motifAnnotations_hgnc <- motifAnnotations

# 4) Initialize SCENIC ----
print('Initialize SCENIC')
scenicOptions <- initializeScenic(
  org = "hgnc",
  dbDir = "/filepath/step15_tfh_subcluster/SCENIC/",
  dbs = c(
    "modified_hg38_10kb_up_down.feather",
    "modified_hg38_500bp_up_100bp_down.feather"
  ),
  nCores = 16
)

# 5) Get RNA counts and metadata from Tfh Seurat object ----
print('Importing Tfh Seurat object')
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

print('Preparing SCENIC inputs from Seurat object')
l3_teaseq_tfh_obj <- JoinLayers(l3_teaseq_tfh_obj, assay = 'RNA') # Required to join RNA layers
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot' # Level 3 analysis Tfh type annotations
exprMat <- l3_teaseq_tfh_obj@assays$RNA$counts
cellInfo <- l3_teaseq_tfh_obj@meta.data['tfh_wnn_annot']
exprMat <- as.matrix(exprMat)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo # Add metadata to SCENIC object

# 6) Filter and select genes to create co-expression network ----
print('Filter and select genes to create co-expression network')
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01) # Default filter
exprMat_filtered <- exprMat[genesKept, ]

# 7) Create co-expression network ----
print('Creating co-expression network')
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 

# 8) Run GENIE3 GRN modeling  ----
print('Running GENIE3 GRN Modeling')
runGenie3(exprMat_filtered_log, scenicOptions)

# GENIE3 step is long - save files at this intermediate step after to avoid restarting from beginning if needed
saveRDS(scenicOptions, '/filepath/step15_tfh_subcluster/SCENIC/scenicOptions_post_Genie3.rds')
saveRDS(exprMat_filtered_log, '/filepath/step15_tfh_subcluster/SCENIC/exprMat_filtered_log_post_Genie3.rds')
saveRDS(exprMat_filtered, '/filepath/step15_tfh_subcluster/SCENIC/exprMat_filtered_post_Genie3.rds')
saveRDS(exprMat, '/filepath/step15_tfh_subcluster/SCENIC/exprMat_post_Genie3.rds')
saveRDS(genesKept, '/filepath/step15_tfh_subcluster/SCENIC/genesKept_post_Genie3.rds')
saveRDS(cellInfo, '/filepath/step15_tfh_subcluster/SCENIC/cellInfo_post_Genie3.rds')

# 9) Build and score GRN model ----

print('Running SCENIC Step 1 - Creating Coexpression Network')
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, '/filepath/step15_tfh_subcluster/SCENIC/scenicOptions1.rds')

print('Running SCENIC Step 2 - Creating Regulons')
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod= NULL) 
saveRDS(scenicOptions, '/filepath/step15_tfh_subcluster/SCENIC/scenicOptions2.rds')

print('Running SCENIC Step 3 - Score Cells')
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, '/filepath/step15_tfh_subcluster/SCENIC/scenicOptions3.rds')

# 10) Binarize regulon activity (optional, used Step 3 AUC scores throughout analyses) ----

print('Running SCENIC Step 4 - Binarize Regulons')
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, '/filepath/step15_tfh_subcluster/SCENIC/scenicOptions4.rds')

print('SCENIC R script finished')