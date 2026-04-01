# Xenium Prime Data Preprocessing Step 8B
# Project dimensionality reduction and subclustering results from downsampled 'sketch' assay of L2 nnCD4 T cells back to full 'RNA' dataset
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

print('Starting R script')

# Set up R working environment

# Set seed
set.seed(26)

# Set working directory
setwd("/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/")

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(BPCells) # 0.3.0
library(presto) # 1.0.0

# Increase max object size allowed for parallelization
options(future.globals.maxSize = 200 * 1024^3)

# Import L2 nnCD4 T cell object from Step 8A Sketch analysis 
print('Importing object')
xp_nncd4t_obj <- readRDS(file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/l2_xp_nncd4_obj_step8a.rds')

# Project L2 Sketch analysis results back to full dataset

# Create reference metadata list of L2 cluster resolutions
sketch_refdata_list <- list(
  full_nncd4t_clust_0.8_res = "sketch_nncd4t_clust_0.8_res",
  full_nncd4t_clust_0.9_res = "sketch_nncd4t_clust_0.9_res",
  full_nncd4t_clust_1.0_res = "sketch_nncd4t_clust_1.0_res"
)

print('Projecting Sketch analysis results back to full dataset')
xp_nncd4t_obj <- ProjectData(
  object = xp_nncd4t_obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:15,
  refdata = sketch_refdata_list
)
print('Sketch projected back to full dataset')

# Save L2 object with sketch results projected back to full dataset
print('Saving object')
saveRDS(xp_nncd4t_obj, file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/l2_xp_nncd4_obj_step8b.rds')
print('Projected sketch object saved')

print('R script done')