# Xenium Prime Data Preprocessing Step 4
# Create Sketch downsampled assay in merged object, regress out probe counts and unique feature counts, compute 30PC UMAP embedding, cluster at several resolutions, and find cluster DEG at each level of resolution
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

# Set working directory
setwd("/filepath/xenium_data_processing/step4_sketch/")

# Increase max object size allowed for parallelization 
options(future.globals.maxSize = 200 * 1024^3)

# Import Seurat object of all six merged samples
print('importing object')
xp.merged.obj <- readRDS('/filepath/xenium_data_processing/step3_merge/xp_obj_step3.rds')

# Create Sketch of full dataset
# SCTransform not recommended for Sketch workflow by Satija Lab - (https://github.com/satijalab/seurat/issues/7336#issuecomment-1555121689)
# Instead, we perform standard NormalizeData > FindVariableFeatures workflow then SketchData
DefaultAssay(xp.merged.obj) <- 'RNA'
print('Normalizing RNA counts')
xp.merged.obj <- NormalizeData(xp.merged.obj, assay = 'RNA')
print('Finding variable RNA features')
xp.merged.obj <- FindVariableFeatures(xp.merged.obj, assay = 'RNA')
print('Calculating Sketch leverage scores')
xp.merged.obj <- SketchData(
  object = xp.merged.obj,
  ncells = 10000, # 10e3 cells per donor, 60e3 total
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# Sketch assay analysis
print('Regression, PCA dimensionality reduction, and UMAP embedding with 30 PC')
DefaultAssay(xp.merged.obj) <- "sketch"
xp.merged.obj <- FindVariableFeatures(xp.merged.obj, assay = 'sketch')
xp.merged.obj <- ScaleData(xp.merged.obj, assay = 'sketch', vars.to.regress = c('nCount_RNA','nFeature_RNA'))
xp.merged.obj <- RunPCA(xp.merged.obj, npcs = 50, reduction.name = 'pca.sketch', reduction.key = 'sketchPC_', future.seed = TRUE)
xp.merged.obj <- FindNeighbors(xp.merged.obj, assay = 'sketch', reduction = "pca.sketch", dims = 1:30, graph.name = c('sketch_nn','sketch_snn'))
xp.merged.obj <- RunUMAP(xp.merged.obj, dims = 1:30, assay = 'sketch', reduction = "pca.sketch", reduction.name = "umap.sketch", reduction.key = "sketchUMAP_", return.model = TRUE)

# Clustering across several levels of resolution
print('Clustering across several levels of resolution')
DefaultAssay(xp.merged.obj) <- 'sketch'
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.15, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.15_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.15_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.2, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.2_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.2_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.3, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.3_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.3_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.35, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.35_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.35_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.4, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.4_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.4_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.5, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.5_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.5_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.6, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.6_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.6_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.7, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.7_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.7_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.8, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.8_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.8_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 0.9, assay = 'sketch', cluster.name = "sketch_bulk_clust_0.9_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_0.9_res <- Idents(xp.merged.obj)
xp.merged.obj <- FindClusters(xp.merged.obj, resolution = 1.0, assay = 'sketch', cluster.name = "sketch_bulk_clust_1.0_res", graph.name = 'sketch_snn')
xp.merged.obj$sketch_bulk_clust_1.0_res <- Idents(xp.merged.obj)

# Save Sketch object with several resolutions of clustering
print('Saving Sketch object with several resolutions of clustering')
saveRDS(xp.merged.obj, file = '/filepath/xenium_data_processing/step4_sketch/xp_obj_step4.rds')
print('Sketch object saved')

# Join RNA layers for cluster annotation
print('Joining RNA layers for finding cluster DEG across resolutions')
DefaultAssay(xp.merged.obj) <- 'sketch'
xp.merged.obj <- JoinLayers(xp.merged.obj, assay = 'sketch')

# Find markers for 0.15 cluster resolution
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.15_res'
rna_markers_res_0.15_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.15_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.15_sketch_regr.rds')

# Find markers for 0.2 cluster resolution
print('Finding cluster markers for resolution 0.2')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.2_res'
rna_markers_res_0.2_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.2_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.2_sketch_regr.rds')

# Find markers for 0.3 cluster resolution
print('Finding cluster markers for resolution 0.3')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.3_res'
rna_markers_res_0.3_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.3_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.3_sketch_regr.rds')

# Find markers for 0.35 cluster resolution
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.35_res'
rna_markers_res_0.35_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.35_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.35_sketch_regr.rds')

# Find markers for 0.4 cluster resolution
print('Finding cluster markers for resolution 0.4')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.4_res'
rna_markers_res_0.4_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.4_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.4_sketch_regr.rds')

# Find markers for 0.5 cluster resolution
print('Finding cluster markers for resolution 0.5')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.5_res'
rna_markers_res_0.5_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.5_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.5_sketch_regr.rds')

# Find markers for 0.6 cluster resolution
print('Finding cluster markers for resolution 0.6')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.6_res'
rna_markers_res_0.6_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.6_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.6_sketch_regr.rds')

# Find markers for 0.7 cluster resolution
print('Finding cluster markers for resolution 0.7')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.7_res'
rna_markers_res_0.7_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.7_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.7_sketch_regr.rds')

# Find markers for 0.8 cluster resolution
print('Finding cluster markers for resolution 0.8')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.8_res'
rna_markers_res_0.8_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.8_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.8_sketch_regr.rds')

# Find markers for 0.9 cluster resolution
print('Finding cluster markers for resolution 0.9')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_0.9_res'
rna_markers_res_0.9_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_0.9_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.9_sketch_regr.rds')

# Find markers for 1.0 cluster resolution
print('Finding cluster markers for resolution 1.0')
Idents(xp.merged.obj) <- 'sketch_bulk_clust_1.0_res'
rna_markers_res_1.0_sketch <- FindAllMarkers(xp.merged.obj, assay = 'sketch')
saveRDS(rna_markers_res_1.0_sketch, '/filepath/xenium_data_processing/step4_sketch/rna_markers_res_1.0_sketch_regr.rds')

print('R script done')