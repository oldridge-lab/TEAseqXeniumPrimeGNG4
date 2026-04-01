# Xenium Prime Data Preprocessing Step 2
# Quality control filtering for each sample (in addition to filtering already performed by Xenium Onboard Analysis)

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC657B, 8YO F with SDB (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Set working directory
setwd("/filepath/xenium_data_processing/step2_qc/")

# Import preprocessed Seurat object
tc657b.obj <- readRDS('/filepath/xenium_data_processing/step1_preprocess/tc657b.preproc.seurat.obj.rds')

# Explore probe count and feature distribution
summary(tc657b.obj$nCount_RNA)
quantile(tc657b.obj$nCount_RNA,0.005) # 25
quantile(tc657b.obj$nCount_RNA,0.01) # 41
quantile(tc657b.obj$nCount_RNA,0.99) # 1357
quantile(tc657b.obj$nCount_RNA,0.995) # 1549.27
summary(tc657b.obj$nFeature_RNA)
quantile(tc657b.obj$nFeature_RNA,0.005) # 23
quantile(tc657b.obj$nFeature_RNA,0.01) # 38
quantile(tc657b.obj$nFeature_RNA,0.99) # 794
quantile(tc657b.obj$nFeature_RNA,0.995) # 865

# Plot unique panel probe features versus total probe counts
FeatureScatter(tc657b.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', log = TRUE) + geom_vline(xintercept = 50) + geom_hline(yintercept = 45) + NoLegend() + labs(title = 'TC657B Pre-QC')

# Determine cells filtered by unique feature and total probe count metrics
sum(tc657b.obj$nCount_RNA < 50) # 16798 cells filtered
sum(tc657b.obj$nCount_RNA < 50) / length(tc657b.obj$nCount_RNA) * 100 # 1.317752 % cells filtered
sum(tc657b.obj$nFeature_RNA < 45) # 16268 cells filtered
sum(tc657b.obj$nFeature_RNA < 45) / length(tc657b.obj$nFeature_RNA) * 100 # 1.276175 % cells filtered
sum(tc657b.obj$nFeature_RNA < 45 | tc657b.obj$nCount_RNA < 50) # 17199 cells filtered
sum(tc657b.obj$nFeature_RNA < 45 | tc657b.obj$nCount_RNA < 50) / ncol(tc657b.obj) * 100 # 1.349209 % cells filtered

# Filter object based on the total of unique RNA features and probe counts per segmented cell
tc657b.qc <- subset(tc657b.obj, subset = nCount_RNA >= 50 & nFeature_RNA >= 45)

# Save object and metadata after QC filtering
saveRDS(tc657b.qc, file = '/filepath/xenium_data_processing/step2_qc/tc657b_post_qc_obj.rds')