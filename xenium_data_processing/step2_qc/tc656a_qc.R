# Xenium Prime Data Preprocessing Step 2
# Quality control filtering for each sample (in addition to filtering already performed by Xenium Onboard Analysis)

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC656A, 4YO M with OSA (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Set working directory
setwd("/filepath/xenium_data_processing/step2_qc/")

# Import preprocessed Seurat object
tc656a.obj <- readRDS('/filepath/xenium_data_processing/step1_preprocess/tc656a.preproc.seurat.obj.rds')

# Explore probe count and feature distribution
summary(tc656a.obj$nCount_RNA)
quantile(tc656a.obj$nCount_RNA,0.005) # 35
quantile(tc656a.obj$nCount_RNA,0.01) # 53
quantile(tc656a.obj$nCount_RNA,0.99) # 1561 
quantile(tc656a.obj$nCount_RNA,0.995) # 1761 
summary(tc656a.obj$nFeature_RNA)
quantile(tc656a.obj$nFeature_RNA,0.005) # 32 
quantile(tc656a.obj$nFeature_RNA,0.01) # 48 
quantile(tc656a.obj$nFeature_RNA,0.99) # 860 
quantile(tc656a.obj$nFeature_RNA,0.995) # 926 

# Plot unique panel probe features versus total probe counts
FeatureScatter(tc656a.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', log = TRUE) + geom_vline(xintercept = 50) + geom_hline(yintercept = 45) + NoLegend() + labs(title = 'TC656A Pre-QC')

# Determine cells filtered by unique feature and total probe count metrics
sum(tc656a.obj$nCount_RNA < 50) # 7716 cells filtered
sum(tc656a.obj$nCount_RNA < 50) / length(tc656a.obj$nCount_RNA) * 100 # 0.8838478 % cells filtered
sum(tc656a.obj$nFeature_RNA < 45) # 7507 cells filtered
sum(tc656a.obj$nFeature_RNA < 45) / length(tc656a.obj$nFeature_RNA) * 100 # 0.8599074 % cells filtered
sum(tc656a.obj$nFeature_RNA < 45 | tc656a.obj$nCount_RNA < 50) # 7974 cells filtered
sum(tc656a.obj$nFeature_RNA < 45 | tc656a.obj$nCount_RNA < 50) / ncol(tc656a.obj) * 100 # 0.913401 % cells filtered

# Filter object based on the total of unique RNA features and probe counts per segmented cell
tc656a.qc <- subset(tc656a.obj, subset = nCount_RNA >= 50 & nFeature_RNA >= 45)

# Save object and metadata after QC filtering
saveRDS(tc656a.qc, file = '/filepath/xenium_data_processing/step2_qc/tc656a_post_qc_obj.rds')