# Xenium Prime Data Preprocessing Step 2
# Quality control filtering for each sample (in addition to filtering already performed by Xenium Onboard Analysis)

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC654B, 8YO F with SDB (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Set working directory
setwd("/filepath/xenium_data_processing/step2_qc/")

# Import preprocessed Seurat object
tc654b.obj <- readRDS('/filepath/xenium_data_processing/step1_preprocess/tc654b.preproc.seurat.obj.rds')

# Explore probe count and feature distribution
summary(tc654b.obj$nCount_RNA)
quantile(tc654b.obj$nCount_RNA,0.005) # 15 
quantile(tc654b.obj$nCount_RNA,0.01) # 27
quantile(tc654b.obj$nCount_RNA,0.99) # 1137 
quantile(tc654b.obj$nCount_RNA,0.995) # 1325
summary(tc654b.obj$nFeature_RNA)
quantile(tc654b.obj$nFeature_RNA,0.005) # 14
quantile(tc654b.obj$nFeature_RNA,0.01) # 25
quantile(tc654b.obj$nFeature_RNA,0.99) # 711.79
quantile(tc654b.obj$nFeature_RNA,0.995) # 789

# Plot unique panel probe features versus total probe counts
FeatureScatter(tc654b.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', log = TRUE) + geom_vline(xintercept = 50) + geom_hline(yintercept = 45) + NoLegend() + labs(title = 'TC654B Pre-QC')

# Determine cells filtered by unique feature and total probe count metrics
sum(tc654b.obj$nCount_RNA < 50) # 20864 cells filtered
sum(tc654b.obj$nCount_RNA < 50) / length(tc654b.obj$nCount_RNA) * 100 # 2.078249 % cells filtered
sum(tc654b.obj$nFeature_RNA < 45) # 20214 cells filtered
sum(tc654b.obj$nFeature_RNA < 45) / length(tc654b.obj$nFeature_RNA) * 100 # 2.013503 % cells filtered
sum(tc654b.obj$nFeature_RNA < 45 | tc654b.obj$nCount_RNA < 50) # 21273 cells filtered
sum(tc654b.obj$nFeature_RNA < 45 | tc654b.obj$nCount_RNA < 50) / ncol(tc654b.obj) * 100 # 2.118989 % cells filtered

# Filter object based on the total of unique RNA features and probe counts per segmented cell
tc654b.qc <- subset(tc654b.obj, subset = nCount_RNA >= 50 & nFeature_RNA >= 45)

# Save object and metadata after QC filtering
saveRDS(tc654b.qc, file = '/filepath/xenium_data_processing/step2_qc/tc654b_post_qc_obj.rds')