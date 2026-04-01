# Step 6 - Remove Empty Droplets (GEM8)

# Initial permissive quality control to remove low UMI count droplets across modalities before performing any HTO demultiplexing or multiplet filtering

# Substeps ----

# 1) Set up R working environment
# 2) Import preprocessed Seurat object for GEM, including RNA, ATAC, ADT, HTO, and SNP modalities.
# 3) Calculate empty droplet filtering thresholds for RNA, ATAC, ADT, and HTO modalities
# 4) Visualize filtering thresholds
# 5) Filter out empty droplets
# 6) Save filtered Seurat object
# 7) Repeat for all 8 GEMs and proceed to next step

# 1) Set up R working environment ----

# Set seed
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

# Set working directory to processedData subfolder for this step
setwd("/filepath/step6_remove_empty")

# 2) Import pre-processed Seurat object for GEM well ----
gem8.preproc <- readRDS('/filepath/step5_preprocess_gems/gem8_preprocessed.rds')

# 3) Calculate empty droplet filtering thresholds for RNA, ATAC, ADT, and HTO modalities ----

# RNA - UMI and Feature thresholds
quantile(gem8.preproc$nCount_RNA, 0.02)
rna.umi.thresh <- 350
quantile(gem8.preproc$nFeature_RNA, 0.02)
rna.feat.thresh <- 250

# ATAC - UMI and Feature thresholds
quantile(gem8.preproc$nCount_ATAC, 0.02)
atac.umi.thresh <- 700
quantile(gem8.preproc$nFeature_ATAC, 0.02)
atac.feat.thresh <- 350

# ADT - only UMI threshold
quantile(gem8.preproc$nCount_ADT, 0.02)
adt.umi.thresh <- 175

# HTO - only UMI threshold
quantile(gem8.preproc$nCount_HTO, 0.02)
hto.umi.thresh <- 50

# Mitochondrial UMI % threshold
DefaultAssay(gem8.preproc) <- 'RNA'
gem8.preproc[["percent.mt"]] <- PercentageFeatureSet(gem8.preproc, pattern = "^MT-") 
quantile(gem8.preproc$percent.mt, 0.98)
mito.thresh <- 30

# 4) Visualize filtering thresholds ----

VlnPlot(gem8.preproc, features = 'nCount_RNA', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = rna.umi.thresh)
VlnPlot(gem8.preproc, features = 'nFeature_RNA', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = rna.feat.thresh)
VlnPlot(gem8.preproc, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = atac.umi.thresh)
VlnPlot(gem8.preproc, features = 'nFeature_ATAC', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = atac.feat.thresh)
VlnPlot(gem8.preproc, features = 'nCount_ADT', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = adt.umi.thresh)
VlnPlot(gem8.preproc, features = 'nCount_HTO', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = hto.umi.thresh)
VlnPlot(gem8.preproc, features = 'percent.mt', pt.size = 0.1, log = TRUE) + geom_hline(yintercept = mito.thresh)

# 5) Filter out empty droplets ----

gem8.filt <- subset(gem8.preproc, subset = 
                    nCount_RNA > rna.umi.thresh &
                    nFeature_RNA > rna.feat.thresh &
                    nCount_ATAC > atac.umi.thresh &
                    nFeature_ATAC > atac.feat.thresh &
                    nCount_ADT > adt.umi.thresh &
                    nCount_HTO > hto.umi.thresh &
                    percent.mt < mito.thresh)

ncol(gem8.preproc) # 5586 droplets initial
ncol(gem8.preproc) - ncol(gem8.filt) # 427 removed
ncol(gem8.filt) # 5159 cell droplets

# 6) Save filtered Seurat object ----
saveRDS(gem8.filt,'/filepath/step6_remove_empty/gem8.filt.rds')

# 7) Repeat for all 8 GEMs and proceed to next step