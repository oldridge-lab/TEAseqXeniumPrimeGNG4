# Xenium Prime Data Preprocessing Step 7
# Level 1 annotation of cell types resolved after projecting 'sketch' assay information to full dataset of ~6.3 million cells
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(BPCells) # 0.3.0
library(presto) # 1.0.0
library(tidyverse) # 2.0.0

# Set working directory
setwd("/filepath/xenium_data_processing/step7_l1_annotation/")

# Import Seurat object after projecting dimensionality reduction and clustering information from downsampled 'sketch' assay to full 'RNA' assay
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step6_project_sketch/xp_obj_step6.rds')

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

# Level 1 annotations for 26 clusters from FindMarkers resolution 1.0
xp_l1_annot <- c(
  '0'  = 'Naive B', # Relative to MBC, less ITGAX 
  '1'  = 'Naive T', # Mixed naive CD4 and CD8 T cell features
  '2'  = 'Mono/Mac', # Mixed monocyte and mature macrophage-related genes
  '3'  = 'FRC', # Fibroblastic reticular cells
  '4'  = 'nnCD4', # non-naive CD4 T cells, Tfh-like signature including TOX2, BCL6, and GNG4
  '5'  = 'ASC', # Antibody-secreting cell phenotype
  '6'  = 'Int Epi', # Found to be spatially positioned between surface and basal epithelium including within crypts using Xenium Explorer
  '7'  = 'LZ GCB', # relative to DZ-like cluster, more CXCR5, LMO2, FCER2 (CD23), CD83, CD40, CD86, CD72, SEMA7A, BCL2A1, EGR3, CIITA
  '8'  = 'DZ GCB', # relative to LZ-like cluster, more CXCR4, AICDA, HMMR, TCF3, POLH, LIG4, and proliferation-related genes
  '9'  = 'Cytotox', # Including features of NK cells and CD8 T cell
  '10' = 'Surf Epi', # Found to be spatially positioned at periphery of epithelium including within crypts using Xenium Explorer
  '11' = 'Cycling', # Largely proliferation-related genes, but also features of several lineages, especially B and T cells
  '12' = 'Basal Epi', # Found to be spatially positioned at basal layer of epithelium including within crypts using Xenium Explorer
  '13' = 'Gran', # Granulocytes, particularly neutrophil-related genes in contrast to c25 Mast-like cells
  '14' = 'Ven Endo', # Venous-like endothelium
  '15' = 'Mural', # Mural-like cells found to be positioned around or near endothelial cells using Xenium Explorer
  '16' = 'Crypt Epi', # Epithelial cell cluster primarily found in crypts but not exterior epithelium using Xenium Explorer 
  '17' = 'DC', # Dendritic cells. Relative to c24 eTAC-like cells, more CPVL ITGAX AOAH SAMHD1 ITGB2 CLEC10A CLEC4A CYBB FOS AHNAK NAIP ANPEP SCIMP CD4 FCGR2A CLEC7A ANXA1 MYCL CIITA
  '18' = 'PDC', # Plasmacytoid dendritic cells
  '19' = 'MBC', # Memory B cells
  '20' = 'Art Endo', # Arterial-like endothelial cells
  '21' = 'FDC', # Follicular dendritic cells
  '22' = 'LEC', # Lymphatic endothelial-like cells
  '23' = 'Trab Fib', # Fibroblasts enriched in trabecular positioning i.e. thick connective tissue bands running through tonsil samples
  '24' = 'eTAC', # Extrathymic AIRE-expressing antigen-presenting cells. Relative to c17 DC cells, more AIRE AMP3 FSCN1 CD83 RGS1 CDKN1A IL7R ADAM12 LAD1 CCL22 DUSP5 KIF2A NRXN2 LY75.
  '25' = 'Mast' # Separate granulocyte cluster from 'Gran' with mast cell signature
)

# For data and rationale supporting cluster annotations, refer to Fig. S14 and Data File S4 for key DEGs, Fig. S16 for spatial distribution, and Supplementary Materials for references and discussion

# Apply L1 annotations to full object
DefaultAssay(xp.obj) <- 'RNA'
xp.obj$l1_annot <- xp.obj$full_bulk_clust_1.0_res
Idents(xp.obj) <- 'l1_annot'
xp.obj <- RenameIdents(xp.obj, xp_l1_annot)
xp.obj$l1_annot <- Idents(xp.obj)
Idents(xp.obj) <- 'l1_annot'

# Assign colors to all 26 clusters for visualization
l1_xp_annot_colors <- c(
  "Naive B" = "#7b7f9c",
  "Naive T" = "#d699ac",
  "Mono/Mac" = "#c49980",
  "FRC" = "#90a6bf",
  "nnCD4" = "#c9744d", 
  "ASC" = "#b37cbf",
  "Int Epi" = "#e49a78",
  "LZ GCB" = "#f7f1e0", 
  "DZ GCB" = "#657ab0",
  "Cytotox" = "#dfcc78", 
  "Surf Epi" = "#d4b5e4", 
  "Cycling" = "#f4a6c7",
  "Basal Epi" = "#e1d4c7", 
  "Gran" = "#728782",
  "Ven Endo" = "#8ab1cc", 
  "Mural" = "#b6b7ba", 
  "Crypt Epi" = "#94a890",
  "DC" = "#9f5a6a",
  "PDC" = "#bad4f4",
  "MBC" = "#971a86", 
  "Art Endo" = "#d66d97",
  "FDC" = "#add5df",
  "LEC" = "#50664f",
  "Trab Fib" = "#e0c1b4", 
  "eTAC" = "#d19a11",
  "Mast" = '#69c9d1'
)

# Compare 'sketch' versus 'RNA' assay UMAP
Idents(xp.obj) <- 'l1_annot'
DefaultAssay(xp.obj) <- 'RNA'
DimPlot(xp.obj, group.by = 'l1_annot', cols = l1_annot_colors, alpha = 0.3, label = TRUE, label.box = TRUE, label.size = 2.5, repel = FALSE, raster = TRUE) + 
  NoLegend() + coord_fixed() + xlab('UMAP1') + ylab('UMAP2') + labs(title = 'Tonsil Level 1 Cell Types [Xenium Prime]')

Idents(xp.obj) <- 'l1_annot'
DefaultAssay(xp.obj) <- 'sketch'
DimPlot(xp.obj, group.by = 'l1_annot', cols = l1_annot_colors, alpha = 0.5, pt.size = 0.8, label = TRUE, label.box = TRUE, label.size = 3, repel = FALSE, raster = FALSE) + 
  NoLegend() + coord_fixed() + xlab('UMAP1') + ylab('UMAP2') + labs(title = 'Tonsil Level 1 Cell Types [Xenium Prime]')

# Import L1 tonsil cluster markers from Step 4 Sketch analysis (nCount/nFeature regression, 30 PC, 1.0 clustering resolution, Sketch workflow)
rna_markers_res_1.0_sketch <- readRDS('/filepath/xenium_data_processing/step4_sketch/rna_markers_res_1.0_sketch_regr.rds')

# Save L1 Sketch assay cluster markers for Data File S4
Idents(xp.obj) <- 'sketch_bulk_clust_1.0_res'
rna_markers_res_1.0_sketch$l1_xp_cluster <- xp_l1_annot[as.character(rna_markers_res_1.0_sketch$cluster)]
l1_xp_clust_rna_sketch_markers <- rna_markers_res_1.0_sketch %>%
  rename( # dataframe columns renamed for clarity
    cluster_number = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_xp_clust_rna_sketch_markers <- l1_xp_clust_rna_sketch_markers %>% arrange(l1_xp_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_xp_clust_rna_sketch_markers <- l1_xp_clust_rna_sketch_markers %>% relocate(c('l1_xp_cluster','cluster_number','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_xp_clust_rna_sketch_markers, '/filepath/xenium_data_processing/step7_l1_annotation/l1_xp_clust_rna_sketch_markers.rds')

# Determine L1 cluster markers using the 'RNA' assay of all ~6.3 million cells (in contrast to downsampled 603e3 cell 'sketch' assay above)
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'full_bulk_clust_1.0_res' # projected sketch-to-RNA assay clustering results from Step 6
xp.obj <- JoinLayers(xp.obj, assay = 'RNA') # joining assay layers to compute DEG
l1_xp_clust_res1_rna_deg <- FindAllMarkers(xp.obj, assay = 'RNA') # default parameters used
l1_xp_clust_res1_rna_deg$xp_l1_cluster <- xp_l1_annot[as.character(l1_xp_clust_res1_rna_deg$cluster)]
l1_xp_rna_assay_markers <- l1_xp_clust_res1_rna_deg %>%
  rename( # dataframe columns renamed for clarity
    cluster_number = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_xp_rna_assay_markers <- l1_xp_rna_assay_markers %>% arrange(xp_l1_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_xp_rna_assay_markers <- l1_xp_rna_assay_markers %>% relocate(c('xp_l1_cluster','cluster_number','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_xp_rna_assay_markers, "/filepath/xenium_data_processing/step7_l1_annotation/l1_xp_rna_assay_markers.rds")

# Export L1 cluster markers for both 'sketch' and 'RNA' assays as tabs of single spreadsheet
clust_marker_df_dir   <- '/filepath/xenium_data_processing/step7_l1_annotation/xenium_l1_cluster_markers'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/xenium_data_processing/step7_l1_annotation/xenium_l1_cluster_markers.xlsx")

# Save Seurat object with L1 cluster annotations
saveRDS(xp.obj, '/filepath/xenium_data_processing/step7_l1_annotation/xp_obj_step7.rds')