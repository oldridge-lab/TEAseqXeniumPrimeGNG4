# Xenium Prime Data Preprocessing Step 2
# Quality control filtering for each sample (in addition to filtering already performed by Xenium Onboard Analysis)

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC659A, 4YO M with SDB (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Set working directory
setwd("/filepath/xenium_data_processing/step2_qc/")

# Import preprocessed Seurat object
tc659a.obj <- readRDS('/filepath/xenium_data_processing/step1_preprocess/tc659a.preproc.seurat.obj.rds')

# Explore probe count and feature distribution
summary(tc659a.obj$nCount_RNA)
quantile(tc659a.obj$nCount_RNA,0.005) # 14
quantile(tc659a.obj$nCount_RNA,0.01) # 26
quantile(tc659a.obj$nCount_RNA,0.99) # 1219 
quantile(tc659a.obj$nCount_RNA,0.995) # 1427.4 
summary(tc659a.obj$nFeature_RNA)
quantile(tc659a.obj$nFeature_RNA,0.005) # 13
quantile(tc659a.obj$nFeature_RNA,0.01) # 24
quantile(tc659a.obj$nFeature_RNA,0.99) # 744 
quantile(tc659a.obj$nFeature_RNA,0.995) # 825 

# Plot unique panel probe features versus total probe counts
FeatureScatter(tc659a.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', log = TRUE) + geom_vline(xintercept = 50) + geom_hline(yintercept = 45) + NoLegend() + labs(title = 'TC659A Pre-QC')

# Determine cells filtered by unique feature and total probe count metrics
sum(tc659a.obj$nCount_RNA < 50) # 28430 cells filtered
sum(tc659a.obj$nCount_RNA < 50) / length(tc659a.obj$nCount_RNA) * 100 # 2.611801 % cells filtered
sum(tc659a.obj$nFeature_RNA < 45) # 27043 cells filtered
sum(tc659a.obj$nFeature_RNA < 45) / length(tc659a.obj$nFeature_RNA) * 100 # 2.48438 % cells filtered
sum(tc659a.obj$nFeature_RNA < 45 | tc659a.obj$nCount_RNA < 50) # 28925 cells filtered
sum(tc659a.obj$nFeature_RNA < 45 | tc659a.obj$nCount_RNA < 50) / ncol(tc659a.obj) * 100 # 2.657275 % cells filtered

# Filter object based on the total of unique RNA features and probe counts per segmented cell
tc659a.qc <- subset(tc659a.obj, subset = nCount_RNA >= 50 & nFeature_RNA >= 45)

# Save object and metadata after QC filtering
saveRDS(tc659a.qc, file = '/filepath/xenium_data_processing/step2_qc/tc659a_post_qc_obj.rds')