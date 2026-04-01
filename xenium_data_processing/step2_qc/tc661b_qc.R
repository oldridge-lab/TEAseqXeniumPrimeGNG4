# Xenium Prime Data Preprocessing Step 2
# Quality control filtering for each sample (in addition to filtering already performed by Xenium Onboard Analysis)

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC661B, 8YO F with SDB (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Set working directory
setwd("/filepath/xenium_data_processing/step2_qc/")

# Import preprocessed Seurat object
tc661b.obj <- readRDS('/filepath/xenium_data_processing/step1_preprocess/tc661b.preproc.seurat.obj.rds')

# Explore probe count and feature distribution
summary(tc661b.obj$nCount_RNA)
quantile(tc661b.obj$nCount_RNA,0.005) # 21 
quantile(tc661b.obj$nCount_RNA,0.01) # 31
quantile(tc661b.obj$nCount_RNA,0.99) # 1177 
quantile(tc661b.obj$nCount_RNA,0.995) # 1349
summary(tc661b.obj$nFeature_RNA)
quantile(tc661b.obj$nFeature_RNA,0.005) # 20
quantile(tc661b.obj$nFeature_RNA,0.01) # 29
quantile(tc661b.obj$nFeature_RNA,0.99) # 727
quantile(tc661b.obj$nFeature_RNA,0.995) # 797

# Plot unique panel probe features versus total probe counts
FeatureScatter(tc661b.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', log = TRUE) + geom_vline(xintercept = 50) + geom_hline(yintercept = 45) + NoLegend() + labs(title = 'TC661B Pre-QC')

# Determine cells filtered by unique feature and total probe count metrics
sum(tc661b.obj$nCount_RNA < 50) # 27133 cells filtered
sum(tc661b.obj$nCount_RNA < 50) / length(tc661b.obj$nCount_RNA) * 100 # 2.445064 % cells filtered
sum(tc661b.obj$nFeature_RNA < 45) # 25605 cells filtered
sum(tc661b.obj$nFeature_RNA < 45) / length(tc661b.obj$nFeature_RNA) * 100 # 2.30737 % cells filtered
sum(tc661b.obj$nFeature_RNA < 45 | tc661b.obj$nCount_RNA < 50) # 27637 cells filtered
sum(tc661b.obj$nFeature_RNA < 45 | tc661b.obj$nCount_RNA < 50) / ncol(tc661b.obj) * 100 # 2.490482 % cells filtered

# Filter object based on the total of unique RNA features and probe counts per segmented cell
tc661b.qc <- subset(tc661b.obj, subset = nCount_RNA >= 50 & nFeature_RNA >= 45)

# Save object and metadata after QC filtering
saveRDS(tc661b.qc, file = '/filepath/xenium_data_processing/step2_qc/tc661b_post_qc_obj.rds')