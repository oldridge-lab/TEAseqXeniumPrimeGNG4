# Xenium Prime Data Preprocessing Step 8A
# Level 2 sketch downsampling and subclustering of 'nnCD4' T cells from L1 object
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

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
options(future.globals.maxSize = 256 * 1024^3)

# Import L1 Seurat object from Step 7
print('Importing L1 clustered Xenium Prime object')
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step7_l1_cluster_annotation/xp_obj_step7.rds')

# Create L2 object of all non-naive CD4 T cells from the L1 clustering analysis
print('Filtering L1 object for nnCD4 T cells')
Idents(xp.obj) <- 'full_bulk_clust_1.0_res'
xp_nncd4t_obj <- subset(xp.obj, idents = '4') # numerical cluster 4 annotated in Step 7 as 'nnCD4'
rm(xp.obj) # to preserve memory

# Downsample L2 nnCD4 T cell subclustering object by standard Sketch workflow
print('Begin preparation of Sketch downsampled object')
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
print('Normalizing RNA counts')
xp_nncd4t_obj <- NormalizeData(xp_nncd4t_obj, assay = 'RNA')
print('Finding variable RNA features')
xp_nncd4t_obj <- FindVariableFeatures(xp_nncd4t_obj, assay = 'RNA')
print('Calculating Sketch leverage scores')
xp_nncd4t_obj <- SketchData(
  object = xp_nncd4t_obj,
  ncells = 10000, # 10e3 cells per donor, 60e3 total
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# L2 Sketch assay analysis
print('Performing dimensionality reduction and clustering of Sketch assay with regression of RNA Counts/Features and resulting 15 PC')
DefaultAssay(xp_nncd4t_obj) <- "sketch"
xp_nncd4t_obj <- FindVariableFeatures(xp_nncd4t_obj, assay = 'sketch')
xp_nncd4t_obj <- ScaleData(xp_nncd4t_obj, assay = 'sketch', vars.to.regress = c('nCount_RNA','nFeature_RNA'))
xp_nncd4t_obj <- RunPCA(xp_nncd4t_obj, npcs = 50, reduction.name = 'pca.sketch', reduction.key = 'sketchPC_', future.seed = TRUE)
xp_nncd4t_obj <- FindNeighbors(xp_nncd4t_obj, assay = 'sketch', reduction = "pca.sketch", dims = 1:15, graph.name = c('sketch_nn','sketch_snn'))
xp_nncd4t_obj <- RunUMAP(xp_nncd4t_obj, dims = 1:15, assay = 'sketch', reduction = "pca.sketch", reduction.name = "umap.sketch", reduction.key = "sketchUMAP_", return.model = TRUE)

# Clustering across several levels of resolution
print('Clustering across several levels of resolution')
DefaultAssay(xp_nncd4t_obj) <- 'sketch'
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.2, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.2_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.2_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.3, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.3_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.3_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.4, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.4_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.4_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.5, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.5_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.5_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.6, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.6_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.6_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.7, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.7_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.7_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.8, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.8_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.8_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 0.9, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_0.9_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_0.9_res <- Idents(xp_nncd4t_obj)
xp_nncd4t_obj <- FindClusters(xp_nncd4t_obj, resolution = 1.0, assay = 'sketch', cluster.name = "sketch_nncd4t_clust_1.0_res", graph.name = 'sketch_snn')
xp_nncd4t_obj$sketch_nncd4t_clust_1.0_res <- Idents(xp_nncd4t_obj)

# Save L2 Sketch object with several resolutions of clustering
print('Saving Sketch object with several resolutions of clustering')
saveRDS(xp_nncd4t_obj, file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/l2_xp_nncd4_obj_step8a.rds')
print('Sketch object saved')

# Join layers to compute DEGs
print('Joining layers to find cluster DEGs across resolutions')
DefaultAssay(xp_nncd4t_obj) <- 'sketch'
xp_nncd4t_obj <- JoinLayers(xp_nncd4t_obj, assay = 'sketch')

# Find markers for 0.2 cluster resolution
print('Finding cluster markers for resolution 0.2')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.2_res'
rna_markers_res_0.2_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.2_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.2_sketch_regr.rds')

# Find markers for 0.3 cluster resolution
print('Finding cluster markers for resolution 0.3')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.3_res'
rna_markers_res_0.3_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.3_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.3_sketch_regr.rds')

# Find markers for 0.4 cluster resolution
print('Finding cluster markers for resolution 0.4')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.4_res'
rna_markers_res_0.4_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.4_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.4_sketch_regr.rds')

# Find markers for 0.5 cluster resolution
print('Finding cluster markers for resolution 0.5')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.5_res'
rna_markers_res_0.5_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.5_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.5_sketch_regr.rds')

# Find markers for 0.6 cluster resolution
print('Finding cluster markers for resolution 0.6')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.6_res'
rna_markers_res_0.6_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.6_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.6_sketch_regr.rds')

# Find markers for 0.7 cluster resolution
print('Finding cluster markers for resolution 0.7')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.7_res'
rna_markers_res_0.7_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.7_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.7_sketch_regr.rds')

# Find markers for 0.8 cluster resolution
print('Finding cluster markers for resolution 0.8')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.8_res'
rna_markers_res_0.8_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.8_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.8_sketch_regr.rds')

# Find markers for 0.9 cluster resolution
print('Finding cluster markers for resolution 0.9')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_0.9_res'
rna_markers_res_0.9_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_0.9_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.9_sketch_regr.rds')

# Find markers for 1.0 cluster resolution
print('Finding cluster markers for resolution 1.0')
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_1.0_res'
rna_markers_res_1.0_sketch <- FindAllMarkers(xp_nncd4t_obj, assay = 'sketch')
saveRDS(rna_markers_res_1.0_sketch, '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_1.0_sketch_regr.rds')

print('R script done')