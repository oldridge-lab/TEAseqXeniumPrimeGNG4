# Xenium Prime Data Preprocessing Step 6
# Project downsampled 'sketch' assay with dimensionality reduction and clustering information for 60e3 cells back to full 'RNA' assay of ~6.3 million cells 
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

print('Starting R script')

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(BPCells) # 0.3.0
library(future) # 1.33.2
library(presto) # 1.0.0
library(tidyverse) # 2.0.0

# Set working directory
setwd("/filepath/xenium_data_processing/step6_project_sketch/")

# Increase max object size allowed for parallelization 
options(future.globals.maxSize = 200 * 1024^3)

# Import Seurat object after performing Sketch workflow with RNA nCount/nFeature regression during feature scaling)
print('Importing object')
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step4_sketch/xp_obj_step4.rds')

print('Projecting Sketch analysis results back to full dataset')

# Create reference metadata list of cluster resolutions
sketch_refdata_list <- list(
  full_bulk_clust_0.35_res = "sketch_bulk_clust_0.35_res",
  full_bulk_clust_0.6_res = "sketch_bulk_clust_0.6_res",
  full_bulk_clust_1.0_res = "sketch_bulk_clust_1.0_res"
)

# Project Sketch analysis results back to full dataset
xp.obj <- ProjectData(
  object = xp.obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = sketch_refdata_list
)
print('Sketch projected back to full dataset')

# Save object with sketch results projected back to full dataset
print('Saving object')
saveRDS(xp.obj, file = '/filepath/xenium_data_processing/step6_project_sketch/xp_obj_step6.rds')
print('Projected sketch object saved')

print('R script done')