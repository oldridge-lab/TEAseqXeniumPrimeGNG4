# Xenium Prime Data Preprocessing Step 2
# Quality control filtering for each sample (in addition to filtering already performed by Xenium Onboard Analysis)

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC653B, 6YO M with SDB (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Set working directory
setwd("/filepath/xenium_data_processing/step2_qc/")

# Import preprocessed Seurat object
tc653b.obj <- readRDS('/filepath/xenium_data_processing/step1_preprocess/tc653b.preproc.seurat.obj.rds')

# Explore probe count and feature distribution
summary(tc653b.obj$nCount_RNA)
quantile(tc653b.obj$nCount_RNA,0.005) # 10
quantile(tc653b.obj$nCount_RNA,0.01) # 17
quantile(tc653b.obj$nCount_RNA,0.99) # 1069.71 
quantile(tc653b.obj$nCount_RNA,0.995) # 1254 
summary(tc653b.obj$nFeature_RNA)
quantile(tc653b.obj$nFeature_RNA,0.005) # 9 
quantile(tc653b.obj$nFeature_RNA,0.01) # 16
quantile(tc653b.obj$nFeature_RNA,0.99) # 658 
quantile(tc653b.obj$nFeature_RNA,0.995) # 731

# Plot unique panel probe features versus total probe counts
FeatureScatter(tc653b.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', log = TRUE) + geom_vline(xintercept = 50) + geom_hline(yintercept = 45) + NoLegend() + labs(title = 'TC653B Pre-QC')

# Determine cells filtered by unique feature and total probe count metrics
sum(tc653b.obj$nCount_RNA < 50) # 36015 cells filtered
sum(tc653b.obj$nCount_RNA < 50) / length(tc653b.obj$nCount_RNA) * 100 # 3.334321 % cells filtered
sum(tc653b.obj$nFeature_RNA < 45) # 34831 cells filtered
sum(tc653b.obj$nFeature_RNA < 45) / length(tc653b.obj$nFeature_RNA) * 100 # 3.224704 % cells filtered
sum(tc653b.obj$nFeature_RNA < 45 | tc653b.obj$nCount_RNA < 50) # 36542 cells filtered
sum(tc653b.obj$nFeature_RNA < 45 | tc653b.obj$nCount_RNA < 50) / ncol(tc653b.obj) * 100 # 3.383111 % cells filtered

# Filter object based on the total of unique RNA features and probe counts per segmented cell
tc653b.qc <- subset(tc653b.obj, subset = nCount_RNA >= 50 & nFeature_RNA >= 45)

# Save object and metadata after QC filtering
saveRDS(tc653b.qc, file = '/filepath/xenium_data_processing/step2_qc/tc653b_post_qc_obj.rds')