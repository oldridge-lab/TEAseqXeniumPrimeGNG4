# Xenium Prime Data Preprocessing Step 8C
# Annotation of L2 nnCD4 T cell subclusters
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

# Set up R working environment ----

# Set seed
set.seed(26)

# Set working directory
setwd("/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/")

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(BPCells) # 0.3.0
library(presto) # 1.0.0
library(tidyverse) # 2.0.0

# Increase max object size allowed for parallelization
options(future.globals.maxSize = 256 * 1024^3)

# Assign colors to numerical clusters before annotation
xp_cols_numb <- c(
  "0"  = "#8ab1cc",
  "1"  = "#add5df",
  "2"  = "#94a890",
  "3"  = "#e0c1b4",
  "4"  = "#e49a78",
  "5"  = "#657ab0",
  "6"  = "#e17794",
  "7"  = "#72d1b4",
  "8"  = "#c9744d",
  "9"  = "#d66d97",
  "10" = "#50664f",
  "11" = "#d4b5e4",
  "12" = "#e1d4c7",
  "13" = "#c9a3c1",
  "14" = "#f4a6c7",
  "15" = "#bad4f4",
  "16" = "#dfcc78",
  "17" = "#728782",
  "18" = "#c49980",
  "19" = "#b6b7ba",
  "20" = "#d699ac",
  "21" = "#50664f",
  "22" = "#4a3c5d",
  "23" = "#72d1b4",
  "24" = "#f7f1e0",
  "25" = "#c9744d",
  "26" = "#b37cbf"
)

# Import Seurat objects ----

# Annotated L1 Seurat object from Step 7
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step7_l1_annotation/xp_obj_step7.rds')

# L2 nnCD4 T cell object from Step 8B (projection of sketch analysis results to full dataset)
xp_nncd4t_obj <- readRDS(file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/l2_xp_nncd4_obj_step8b.rds')

# Inspect L2 nnCD4 T cell object cluster distribution and DEGs at each level of clustering resolution ----

# Inspect number of output clusters for each level of clustering resolution
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.2_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.3_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.4_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.5_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.6_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.7_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.8_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_0.9_res)
unique(xp_nncd4t_obj$sketch_nncd4t_clust_1.0_res) # 16 clusters - ultimately used this resolution for downstream analyses

# Inspect UMAP distribution of clusters at each level of resolution 
DefaultAssay(xp_nncd4t_obj) <- 'sketch'
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.2_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.3_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.4_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.5_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.6_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.7_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.8_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_0.9_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()
DimPlot(xp_nncd4t_obj, reduction = 'umap.sketch', group.by = 'sketch_nncd4t_clust_1.0_res', cols = xp_cols_numb, label = TRUE, label.box = TRUE) + coord_fixed()

# Assess enrichment of genes of interest across clusters
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'CXCR5', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'PDCD1', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'BCL6', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'GNG4', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'IL21', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'POU2AF1', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'NFATC1', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'FOXP3', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'IL10', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'TBX21', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'GATA3', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'RORC', order = TRUE) + coord_fixed()
FeaturePlot(xp_nncd4t_obj, reduction = 'umap.sketch', cols = c('grey','darkred'), features = 'PRDM1', order = TRUE) + coord_fixed()

# Inspect DEGs between clusters at each level of resolution from Sketch analysis in Step 8A
nncd4t_clust_markers_res_0.1_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.1_sketch_regr.rds')
nncd4t_clust_markers_res_0.2_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.2_sketch_regr.rds')
nncd4t_clust_markers_res_0.3_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.3_sketch_regr.rds')
nncd4t_clust_markers_res_0.4_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.4_sketch_regr.rds')
nncd4t_clust_markers_res_0.5_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.5_sketch_regr.rds')
nncd4t_clust_markers_res_0.6_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.6_sketch_regr.rds')
nncd4t_clust_markers_res_0.7_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.7_sketch_regr.rds')
nncd4t_clust_markers_res_0.8_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.8_sketch_regr.rds')
nncd4t_clust_markers_res_0.9_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_0.9_sketch_regr.rds')
nncd4t_clust_markers_res_1.0_sketch <- readRDS('/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/rna_markers_res_1.0_sketch_regr.rds')

# XP L2 nnCD4T Cluster Markers - using full 'RNA' assay dataset rather than downsampled 'sketch' assay
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_nncd4t_obj <- JoinLayers(xp_nncd4t_obj, assay = 'RNA')
Idents(xp_nncd4t_obj) <- 'full_nncd4t_clust_1.0_res' # 16 clusters - ultimately used this resolution for downstream analyses
l2_xp_nncd4t_clust_res1_rna_deg <- FindAllMarkers(xp_nncd4t_obj, assay = 'RNA')

# Annotate L2 nnCD4 T cell subclusters ----

# Refer to Fig. 4J, Fig. S15, and Data File S4 for key DEGs, Fig. S16 for spatial distribution, and Supplementary Materials for references and additional discussion of rationale behind cluster annotations.

l2_nncd4t_annot <- c(
  '0'  = 'Tfh DZ',
  '1'  = 'Tfh Mant',
  '2'  = 'Tcm 1',
  '3'  = 'Tfh CXCL13',
  '4'  = 'Tfh TOX2',
  '5'  = 'Tfh S1PR2',
  '6'  = 'Treg/fr',
  '7'  = 'Tconv Myl', 
  '8'  = 'Tfh CCR6',
  '9'  = 'Trm', 
  '10' = 'Tconv', 
  '11' = 'Tcm 2',
  '12' = 'Tfh Myl',
  '13' = 'Tfh NFATC1',
  '14' = 'Tfh Circ', 
  '15' = 'Tfh PRDM1'
)

# Transfer nnCD4 T cell subcluster annotations to L2 Seurat object
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_nncd4t_obj$l2_nncd4t_annot <- xp_nncd4t_obj$full_nncd4t_clust_1.0_res
Idents(xp_nncd4t_obj) <- 'l2_nncd4t_annot'
xp_nncd4t_obj <- RenameIdents(xp_nncd4t_obj, l2_nncd4t_annot)
xp_nncd4t_obj$l2_nncd4t_annot <- Idents(xp_nncd4t_obj)

# Create grouped labels for Tfh versus nonTfh subclusters
xp_nncd4t_obj_md <- xp_nncd4t_obj@meta.data
xp_nncd4t_obj_md <- xp_nncd4t_obj_md %>%
  mutate(l2_tfh_vs_nontfh = if_else(
    !is.na(l2_nncd4t_annot) & startsWith(as.character(l2_nncd4t_annot), "Tfh"),
    "Tfh", "nonTfh"
  ))
xp_nncd4t_obj@meta.data <- xp_nncd4t_obj_md 
table(xp_nncd4t_obj$l2_tfh_vs_nontfh) # 330,994 cells total - 235811 Tfh (71.2%) versus 95183 nonTfh (28.8%)

# Export nnCD4 T cell subcluster markers for Data File S4 ----

# L2 markers from sketch assay
DefaultAssay(xp_nncd4t_obj) <- 'sketch'
Idents(xp_nncd4t_obj) <- 'sketch_nncd4t_clust_1.0_res' # for 16 clusters
nncd4t_clust_markers_res_1.0_sketch$l2_nncd4t_cluster <- l2_nncd4t_annot[as.character(nncd4t_clust_markers_res_1.0_sketch$cluster)]
l2_nncd4t_xp_rna_sketch_markers <- nncd4t_clust_markers_res_1.0_sketch %>%
  rename(
    cluster_number = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  ) 
l2_nncd4t_xp_rna_sketch_markers <- l2_nncd4t_xp_rna_sketch_markers %>% arrange(l2_nncd4t_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_nncd4t_xp_rna_sketch_markers <- l2_nncd4t_xp_rna_sketch_markers %>% relocate(c('l2_nncd4t_cluster','cluster_number','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_nncd4t_xp_rna_sketch_markers, '/filepath/fig4/fig4_supp/xenium_l2_nncd4t_cluster_markers/l2_nncd4t_xp_rna_sketch_markers.rds')

# L2 markers from full RNA assay
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
Idents(xp_nncd4t_obj) <- 'full_nncd4t_clust_1.0_res' # projected Sketch clustering to full object
l2_xp_nncd4t_clust_res1_rna_deg$l2_nncd4t_cluster <- l2_nncd4t_annot[as.character(l2_xp_nncd4t_clust_res1_rna_deg$cluster)]
l2_nncd4t_xp_rna_assay_markers <- l2_xp_nncd4t_clust_res1_rna_deg %>%
  rename(
    cluster_number = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_nncd4t_xp_rna_assay_markers <- l2_nncd4t_xp_rna_assay_markers %>% arrange(l2_nncd4t_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_nncd4t_xp_rna_assay_markers <- l2_nncd4t_xp_rna_assay_markers %>% relocate(c('l2_nncd4t_cluster','cluster_number','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_nncd4t_xp_rna_assay_markers, '/filepath/fig4/fig4_supp/xenium_l2_nncd4t_cluster_markers/l2_nncd4t_xp_rna_assay_markers.rds')

# Export L2 nnCD4 T cell markers from 'sketch' and 'RNA' assays as single .xlsx spreadsheet
clust_marker_df_dir   <- '/filepath/fig4/fig4_supp/xenium_l2_nncd4t_cluster_markers'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/fig4/fig4_supp/xenium_l2_nncd4t_cluster_markers/xenium_l2_nncd4t_cluster_markers.xlsx")  

# Additional comparisons for cluster annotation - cTfh cluster vs nonTfh group, and Tfh group vs nonTfh group - shown in Data File S4 ----

# Get names of clusters from Tfh vs nonTfh groups
tfh_clust_names <- l2_nncd4t_annot[startsWith(l2_nncd4t_annot, "Tfh")]
nontfh_clust_names <- l2_nncd4t_annot[!startsWith(l2_nncd4t_annot, "Tfh")]

# Join layers for both assays to compute DEGs
xp_nncd4t_obj <- JoinLayers(xp_nncd4t_obj, assay = 'sketch')
xp_nncd4t_obj <- JoinLayers(xp_nncd4t_obj, assay = 'RNA')

# Tfh group vs nonTfh group - Sketch assay
l2_xp_tfh_vs_nontfh_deg_sketch <- FindMarkers(xp_nncd4t_obj, ident.1 = tfh_clust_names, ident.2 = nontfh_clust_names,  assay = 'sketch', min.pct = 0, logfc.threshold = 0)
l2_xp_tfh_vs_nontfh_deg_sketch$gene <- rownames(l2_xp_tfh_vs_nontfh_deg_sketch)
l2_xp_tfh_vs_nontfh_deg_sketch_markers <- l2_xp_tfh_vs_nontfh_deg_sketch %>%
  rename(
    gene_symbol = gene,
    avg_log2fc_tfh_vs_nontfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_ctfh = pct.1,
    pct_pos_in_nontfh = pct.2
  ) 
l2_xp_tfh_vs_nontfh_deg_sketch_markers <- l2_xp_tfh_vs_nontfh_deg_sketch_markers %>% arrange(desc(avg_log2fc_tfh_vs_nontfh), p_val_adj)
l2_xp_tfh_vs_nontfh_deg_sketch_markers <- l2_xp_tfh_vs_nontfh_deg_sketch_markers %>% relocate(c('gene_symbol','avg_log2fc_tfh_vs_nontfh','p_val_raw','p_val_adj'))
saveRDS(l2_xp_tfh_vs_nontfh_deg_sketch_markers, '/filepath/fig4/fig4_supp/xenium_l2_tfh_clust_comparisons/l2_xp_tfh_vs_nontfh_deg_sketch_markers.rds')

# Tfh group vs nonTfh group - RNA assay
l2_xp_tfh_vs_nontfh_deg_rna <- FindMarkers(xp_nncd4t_obj, ident.1 = tfh_clust_names, ident.2 = nontfh_clust_names,  assay = 'RNA', min.pct = 0, logfc.threshold = 0)
l2_xp_tfh_vs_nontfh_deg_rna$gene <- rownames(l2_xp_tfh_vs_nontfh_deg_rna)
l2_xp_tfh_vs_nontfh_deg_rna_markers <- l2_xp_tfh_vs_nontfh_deg_rna %>%
  rename(
    gene_symbol = gene,
    avg_log2fc_tfh_vs_nontfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_ctfh = pct.1,
    pct_pos_in_nontfh = pct.2
  ) 
l2_xp_tfh_vs_nontfh_deg_rna_markers <- l2_xp_tfh_vs_nontfh_deg_rna_markers %>% arrange(desc(avg_log2fc_tfh_vs_nontfh), p_val_adj)
l2_xp_tfh_vs_nontfh_deg_rna_markers <- l2_xp_tfh_vs_nontfh_deg_rna_markers %>% relocate(c('gene_symbol','avg_log2fc_tfh_vs_nontfh','p_val_raw','p_val_adj'))
saveRDS(l2_xp_tfh_vs_nontfh_deg_rna_markers, '/filepath/fig4/fig4_supp/xenium_l2_tfh_clust_comparisons/l2_xp_tfh_vs_nontfh_deg_rna_markers.rds')

# Tfh Circ vs nonTfh - Sketch Assay
Idents(xp_nncd4t_obj) <- 'l2_nncd4t_annot'
tfhcirc_vs_nontfh_deg_sketch <- FindMarkers(xp_nncd4t_obj, ident.1 = 'Tfh Circ', ident.2 = nontfh_clust_names, assay = 'sketch', min.pct = 0, logfc.threshold = 0)
tfhcirc_vs_nontfh_deg_sketch$gene <- rownames(tfhcirc_vs_nontfh_deg_sketch) # CXCR5 increased in cTfh - BHLHE40 decreased
tfhcirc_vs_nontfh_deg_sketch_markers <- tfhcirc_vs_nontfh_deg_sketch %>%
  rename(
    gene_symbol = gene,
    avg_log2fc_ctfh_vs_nontfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_ctfh = pct.1,
    pct_pos_in_nontfh = pct.2
  ) 
tfhcirc_vs_nontfh_deg_sketch_markers <- tfhcirc_vs_nontfh_deg_sketch_markers %>% arrange(desc(avg_log2fc_ctfh_vs_nontfh), p_val_adj)
tfhcirc_vs_nontfh_deg_sketch_markers <- tfhcirc_vs_nontfh_deg_sketch_markers %>% relocate(c('gene_symbol','avg_log2fc_ctfh_vs_nontfh','p_val_raw','p_val_adj'))
saveRDS(tfhcirc_vs_nontfh_deg_sketch_markers, '/filepath/fig4/fig4_supp/xenium_l2_tfh_clust_comparisons/l2_xp_tfhcirc_vs_nontfh_deg_sketch_markers.rds')

# Tfh Circ vs nonTfh - RNA Assay
tfhcirc_vs_nontfh_deg_rna <- FindMarkers(xp_nncd4t_obj, ident.1 = 'Tfh Circ', ident.2 = nontfh_clust_names, assay = 'RNA', min.pct = 0, logfc.threshold = 0)
tfhcirc_vs_nontfh_deg_rna$gene <- rownames(tfhcirc_vs_nontfh_deg_rna) # CXCR5, BCL6, PDCD1, IL21, and CXCL13 increased in cTfh - PRDM1 KLF2 BHLHE40 decreased
tfhcirc_vs_nontfh_deg_rna_markers <- tfhcirc_vs_nontfh_deg_rna %>%
  rename(
    gene_symbol = gene,
    avg_log2fc_ctfh_vs_nontfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_ctfh = pct.1,
    pct_pos_in_nontfh = pct.2
  ) 
tfhcirc_vs_nontfh_deg_rna_markers <- tfhcirc_vs_nontfh_deg_rna_markers %>% arrange(desc(avg_log2fc_ctfh_vs_nontfh), p_val_adj)
tfhcirc_vs_nontfh_deg_rna_markers <- tfhcirc_vs_nontfh_deg_rna_markers %>% relocate(c('gene_symbol','avg_log2fc_ctfh_vs_nontfh','p_val_raw','p_val_adj'))
saveRDS(tfhcirc_vs_nontfh_deg_rna_markers, '/filepath/fig4/fig4_supp/xenium_l2_tfh_clust_comparisons/l2_xp_tfhcirc_vs_nontfh_deg_rna_markers.rds')

# Export additional L2 Tfh vs nonTfh comparison .rds marker files using 'sketch' and 'RNA' assays as sheets of single .xlsx spreadsheet
clust_marker_df_dir   <- '/filepath/fig4/fig4_supp/xenium_l2_tfh_clust_comparisons'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/fig4/fig4_supp/xenium_l2_tfh_clust_comparisons/xenium_l2_tfh_cluster_comparison_markers.xlsx")  

# Build combined L1 cluster and L2 nnCD4 subcluster annotation vector for L1 object ----

# Get L2 nnCD4 T cell cluster annotations per cell
l2_nncd4t_annot <- as.character(xp_nncd4t_obj@meta.data$l2_nncd4t_annot)
names(l2_nncd4t_annot) <- rownames(xp_nncd4t_obj@meta.data)

# Get L1 object metadata, barcodes, and L1 cluster annotations
l1_obj_md  <- xp.obj@meta.data
l1_obj_barcodes <- rownames(l1_obj_md)
l1_obj_labs   <- as.character(l1_obj_md$l1_annot)

# Replace 'nnCD4' annotations in L1 object with L2 annotations, while retaining L1 annotations for all 25 other clusters
l1_obj_md$l1_l2_nncd4t_annot <- ifelse(
  test = l1_obj_labs == "nnCD4",
  yes = l2_nncd4t_annot[l1_obj_barcodes],
  no = l1_obj_labs
)

# Replace L1 object metadata with new metadata including combined L1/L2 annotations
xp.obj@meta.data <- l1_obj_md
table(xp.obj@meta.data$l1_l2_nncd4t_annot)

# Extract combined L1/L2 annotations for visualization in Xenium Explorer ----

# TC653B - extract combined L1/L2 annotations and save as .csv to visualize in Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc653b_l1_with_l2_nncd4t_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC656A - extract combined L1/L2 annotations and save as .csv to visualize in Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
cells_tc656a <- colnames(xp.obj)[ xp.obj$donor == "TC656A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc656a)
idents_tc656a <- Idents(xp.obj)[ cells_tc656a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc656a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc656a_l1_with_l2_nncd4t_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC654B - extract combined L1/L2 annotations and save as .csv to visualize in Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
cells_tc654b <- colnames(xp.obj)[ xp.obj$donor == "TC654B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc654b)
idents_tc654b <- Idents(xp.obj)[ cells_tc654b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc654b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc654b_l1_with_l2_nncd4t_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC661B - extract combined L1/L2 annotations and save as .csv to visualize in Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
cells_tc661b <- colnames(xp.obj)[ xp.obj$donor == "TC661B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc661b)
idents_tc661b <- Idents(xp.obj)[ cells_tc661b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc661b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc661b_l1_with_l2_nncd4t_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC659A - extract combined L1/L2 annotations and save as .csv to visualize in Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
cells_tc659a <- colnames(xp.obj)[ xp.obj$donor == "TC659A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc659a)
idents_tc659a <- Idents(xp.obj)[ cells_tc659a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc659a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc659a_l1_with_l2_nncd4t_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC657B - extract combined L1/L2 annotations and save as .csv to visualize in Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
cells_tc657b <- colnames(xp.obj)[ xp.obj$donor == "TC657B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc657b)
idents_tc657b <- Idents(xp.obj)[ cells_tc657b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc657b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc657b_l1_with_l2_nncd4t_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# Build combined L1 cluster and L2 grouped Tfh vs nonTfh annotation vector for L1 object ----

# Get L2 Tfh vs nonTfh annotation vector
tfh_vs_nontfh_annot_vec <- as.character(xp_nncd4t_obj@meta.data$l2_tfh_vs_nontfh)
names(tfh_vs_nontfh_annot_vec) <- rownames(xp_nncd4t_obj@meta.data)

# Build combined L1/L2 annotation vector, wherein L1 'nnCD4' annotations are replaced with L2 Tfh/nonTfh group label, while all other L1 annotations are retained
l1_obj_md$l1_l2_tfh_vs_nontfh <- ifelse(
  test  = l1_obj_labs == "nnCD4",
  yes = tfh_vs_nontfh_annot_vec[l1_obj_barcodes],
  no = l1_obj_labs
)

# Replace L1 object metadata with updated metadata
xp.obj@meta.data <- l1_obj_md
table(xp.obj@meta.data$l1_l2_tfh_vs_nontfh)

# Extract combined L1 and L2 Tfh vs nonTfh annotations for Xenium Explorer visualization ----

# TC653B - combined L1/L2 Tfh vs nonTfh annotations
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_tfh_vs_nontfh'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc653b_l1_with_l2_tfh_vs_nontfh_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC654B - combined L1/L2 Tfh vs nonTfh annotations
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_tfh_vs_nontfh'
cells_tc654b <- colnames(xp.obj)[ xp.obj$donor == "TC654B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc654b)
idents_tc654b <- Idents(xp.obj)[ cells_tc654b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc654b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc654b_l1_with_l2_tfh_vs_nontfh_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC656A - combined L1/L2 Tfh vs nonTfh annotations
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_tfh_vs_nontfh'
cells_tc656a <- colnames(xp.obj)[ xp.obj$donor == "TC656A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc656a)
idents_tc656a <- Idents(xp.obj)[ cells_tc656a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc656a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc656a_l1_with_l2_tfh_vs_nontfh_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC657B - combined L1/L2 Tfh vs nonTfh annotations
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_tfh_vs_nontfh'
cells_tc657b <- colnames(xp.obj)[ xp.obj$donor == "TC657B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc657b)
idents_tc657b <- Idents(xp.obj)[ cells_tc657b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc657b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc657b_l1_with_l2_tfh_vs_nontfh_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC659A - combined L1/L2 Tfh vs nonTfh annotations
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_tfh_vs_nontfh'
cells_tc659a <- colnames(xp.obj)[ xp.obj$donor == "TC659A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc659a)
idents_tc659a <- Idents(xp.obj)[ cells_tc659a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc659a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc659a_l1_with_l2_tfh_vs_nontfh_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC661B - combined L1/L2 Tfh vs nonTfh annotations
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_tfh_vs_nontfh'
cells_tc661b <- colnames(xp.obj)[ xp.obj$donor == "TC661B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc661b)
idents_tc661b <- Idents(xp.obj)[ cells_tc661b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc661b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/tc661b_l1_with_l2_tfh_vs_nontfh_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# Create L3 Tfh-only subset object from L2 object Tfh subclusters ----

DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_tfh_subset_cell_ids <- rownames(xp_nncd4t_obj@meta.data)[
  grepl("^Tfh", xp_nncd4t_obj@meta.data$l2_nncd4t_annot)
  ]
xp_tfh_obj <- subset(xp_nncd4t_obj, cells = xp_tfh_subset_cell_ids) # 235811 Tfh cells total in L3 object

# No further dimensionality reduction or subclustering performed of L3 object - only created to facilitate Tfh DEG and spatial positioning comparisons in downstream analyses

# Save Seurat objects with updated annotations ----

# L1 object with updated metadata including L2 annotations
saveRDS(xp.obj, file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/xp_l1_obj_step8c.rds')

# Annotated L2 nnCD4 T subclustering object 
saveRDS(xp_nncd4t_obj, file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/xp_l2_nncd4t_obj_step8c.rds')

# New L3 object including only Tfh cells from L2 object - no further dimensionality reduction or clustering performed
saveRDS(xp_tfh_obj, file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/xp_l3_tfh_subset_obj_step8c.rds')