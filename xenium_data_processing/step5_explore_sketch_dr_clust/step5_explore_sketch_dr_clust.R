# Xenium Prime Data Preprocessing Step 5
# Exploring clustering resolutions to project from downsampled 'sketch' assay (60e3 cells) to full 'RNA' assay ( ~6.3 million cells)
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

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
setwd("/filepath/xenium_data_processing/step5_explore_sketch_dr_clust/")

# Increase max object size allowed for parallelization 
options(future.globals.maxSize = 200 * 1024^3)

# Import Seurat object after performing Sketch workflow with RNA nCount/nFeature regression during feature scaling)
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step4_sketch/xp_obj_step4.rds')

# Import downsampled 'sketch' assay cluster markers for each level of resolution determined in Step 4 (with RNA nCount/nFeature regression during feature centering and scaling)
rna_markers_res_0.15_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.15_sketch_regr.rds')
rna_markers_res_0.2_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.2_sketch_regr.rds')
rna_markers_res_0.3_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.3_sketch_regr.rds')
rna_markers_res_0.35_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.35_sketch_regr.rds')
rna_markers_res_0.4_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.4_sketch_regr.rds')
rna_markers_res_0.5_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.5_sketch_regr.rds')
rna_markers_res_0.6_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.6_sketch_regr.rds')
rna_markers_res_0.7_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.7_sketch_regr.rds')
rna_markers_res_0.8_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.8_sketch_regr.rds')
rna_markers_res_0.9_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_0.9_sketch_regr.rds')
rna_markers_res_1.0_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_1.0_sketch_regr.rds')

# Compare clustering at each level of resolution
DefaultAssay(xp.obj) <- 'sketch'
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.15_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.2_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.3_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.35_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.4_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.5_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.6_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.7_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.8_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_0.9_res') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'sketch_bulk_clust_1.0_res') + coord_fixed()

# Assess batch effects across donor demographic and biological variables
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'donor') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'asab') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'age') + coord_fixed()
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'indication') + coord_fixed()

# Assess cell clustering driven by quality control and technical variables
DimPlot(xp.obj, reduction = 'umap.sketch', label = TRUE, label.box = TRUE, group.by = 'slide') + coord_fixed()
Idents(xp.obj) <- 'sketch_bulk_clust_1.0_res' # ultimately selected this resolution for Level 1 clustering and annotation in manuscript
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_RNA')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_RNA', min.cutoff = 'q1', max.cutoff = 'q99')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nFeature_RNA')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nFeature_RNA',  min.cutoff = 'q1', max.cutoff = 'q99')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_ControlProbe')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_GenomicControl')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_BlankCodeword')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_ControlCodeword')
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'nCount_DeprecatedCodeword')

# Ultimately we did not observe any substantial batch effects between the two Xenium Prime slides or between tonsil samples, and therefore did not perform any integration (i.e. by Harmony)

# Example DEGs for determining ideal cluster resolution to separate major expected immune and non-immune cell types in tonsil
Idents(xp.obj) <- 'sketch_bulk_clust_1.0_res' # ultimately selected this resolution for Level 1 clustering and annotation in manuscript
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'CD3E', order = TRUE) + coord_fixed() # T cells
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'MS4A1', order = TRUE) + coord_fixed() # B cells
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'AICDA', order = TRUE) + coord_fixed() # GC B cells
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'MKI67', order = TRUE) + coord_fixed() # GC B cells and proliferating cells
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'AIRE', order = TRUE) + coord_fixed() # c24 eTAC
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'ITGAX', order = TRUE) + coord_fixed() # Myeloid and MBCs
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'DSG3', order = TRUE) + coord_fixed() # Epithelial
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'GNG4', order = TRUE) + coord_fixed() # c4 nnCD4
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'IL21', order = TRUE) + coord_fixed() # c4 nnCD4
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'TOX2', order = TRUE) + coord_fixed() # c4 nnCD4
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = 'CHGB', order = TRUE) + coord_fixed() # c4 nnCD4
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = c('IL1B','CSF3R','PTGS2','PLAUR')) # Neutrophils
FeaturePlot(xp.obj, reduction = 'umap.sketch', label = TRUE, features = c('MMP9','CD68','CTSL','ENPP2')) # Macrophages