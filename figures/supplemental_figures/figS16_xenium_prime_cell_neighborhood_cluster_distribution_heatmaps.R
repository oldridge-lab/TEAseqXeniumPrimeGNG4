# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S16
# Fig. S16 – Immune and non-immune cell types form diverse cellular neighborhoods in tonsils with distinct CD4 T cell subset composition.

# Cell neighborhood (CN) analysis used L1 cluster annotations as input (not L2 T cell subclusters), k = 40 neighbors, and n = 10 neighborhoods
# Original method from Schürch, Christian M. et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. Cell, Volume 182, Issue 5, 1341 - 1359.e19
# Code adapted from related GitHub page - https://github.com/nolanlab/NeighborhoodCoordination
# CN generation for tonsil Xenium Prime dataset was performed in collaboration with Yutong Zhu from Oldridge Lab.
# For details regarding CN generation refer to Python notebooks in the Xenium Data Processing Step 9 folder.
# CNs were annotated as detailed in Xenium Data Preprocessing Step 10 script.

# Set up R working environment ----

# Set working directory
setwd('/filepath/fig4/fig4_supp/')

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Matrix) # 1.7-0
library(tidyselect) # 1.2.1
library(Seurat) # 5.1.0
library(SeuratObject) # 5.0.2
library(data.table) # 1.15.4
library(readr) # 2.1.5
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4
library(stringr) # 1.5.1
library(presto) # 1.0.0
library(paletteer) # 1.6.0
library(RColorBrewer) # 1.1-3
library(EnhancedVolcano) # 1.22.0
library(tidyverse) # 2.0.0
library(BPCells) # 0.3.0
library(tibble) # 3.2.1
library(purrr) # 1.0.2
library(circlize) # 0.4.16
library(writexl) # 1.5.40
library(ggnewscale) # 0.5.0
library(readxl) # 1.4.3
library(openxlsx) # 4.2.8
library(scales) # 1.3.0
library(ggVennDiagram) # 1.5.2
library(pheatmap) # 1.0.12

# Import Xenium Prime L1 and L2 Seurat objects ----

# L1 object - primary immune and non-immune lineages, 26 total clusters
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l1_obj_step10.rds')

# L2 object - nnCD4 T cells, 16 total subclusters of L1 nnCD4 cluster
xp_nncd4t_obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l2_nncd4t_obj_step10.rds')

# Fig. S16G - Heatmap for L1 cluster enrichment across CNs (initial CN analysis version) ----

# Create matrix of L1 cluster cell counts per CN
xp_l1l2clust_cn_df <- xp.obj@meta.data %>%
  dplyr::select(cluster = l1_annot,
                neighborhood = cn_annot)
xp_l1l2clust_cn_table <- table(xp_l1l2clust_cn_df$cluster, xp_l1l2clust_cn_df$neighborhood)
xp_l1l2clust_cn_mtx <- as.matrix(xp_l1l2clust_cn_table)

# Compute within-cluster cell count Z-scores across CNs
xp_l1l2clust_cn_mtx_scaled <- t(scale(t(xp_l1l2clust_cn_mtx)))

# Build color palette for Z-scores centered around 0
cn_heatmap_col_breaks <- seq(-max(abs(xp_l1l2clust_cn_mtx_scaled)), max(abs(xp_l1l2clust_cn_mtx_scaled)), length.out = 100)
cn_heatmap_cols  <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(length(cn_heatmap_col_breaks) - 1))

# Specify CN names for heatmap
colnames(xp_l1l2clust_cn_mtx_scaled) <- c('CN0 Fol Border',
                                          'CN1 TCZ Inner',
                                          'CN2 GC LZ',
                                          'CN3 TCZ Outer',
                                          'CN4 Epi Outer',
                                          'CN5 GZ DZ',
                                          'CN6 Mantle',
                                          'CN7 Epi Inner',
                                          'CN8 ASC Rich',
                                          'CN9 Connective') 

# Transpose matrix to plot clusters as columns and CNs as rows on heatmap
xp_l1l2clust_cn_mtx_scaled_t <- t(xp_l1l2clust_cn_mtx_scaled)

# Assemble heatmap
pdf('/filepath/fig4/fig4_supp/cn_orig_l1_clust_distr_heatmap.pdf', 
    height = 5, width = 10,
    family = "sans")
cn_orig_l1_clust_distr_heatmap <- pheatmap(
  xp_l1l2clust_cn_mtx_scaled_t,
  color         = cn_heatmap_cols,
  breaks        = cn_heatmap_col_breaks,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  fontsize      = 13,
  fontsize_row  = 13,
  fontsize_col  = 13,
  angle_col     = "45",
  border_color  = "grey80",
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend        = TRUE
)
cn_orig_l1_clust_distr_heatmap
dev.off()

# Fig. S16H - Heatmap for L1 cluster enrichment across CNs (merged GC CN analysis version) ----

# Create matrix of L1 cluster cell counts per CN
xp_l1l2clust_cn_df <- xp.obj@meta.data %>%
  dplyr::select(cluster = l1_annot,
                neighborhood = cn_gc_closed_annot)
xp_l1l2clust_cn_table <- table(xp_l1l2clust_cn_df$cluster, xp_l1l2clust_cn_df$neighborhood)
xp_l1l2clust_cn_mtx <- as.matrix(xp_l1l2clust_cn_table)

# Transpose matrix to plot clusters as columns and CNs as rows on heatmap
xp_l1l2clust_cn_mtx_scaled <- t(scale(t(xp_l1l2clust_cn_mtx)))

# Build color palette for Z-scores centered around 0
cn_heatmap_col_breaks <- seq(-max(abs(xp_l1l2clust_cn_mtx_scaled)), max(abs(xp_l1l2clust_cn_mtx_scaled)), length.out = 100)
cn_heatmap_cols  <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(length(cn_heatmap_col_breaks) - 1))

# Specify CN names for heatmap
colnames(xp_l1l2clust_cn_mtx_scaled) <- c('CN0 Fol Border',
                                          'CN1 TCZ Inner',
                                          'CN2 GC',
                                          'CN3 TCZ Outer',
                                          'CN4 Epi Outer',
                                          'CN6 Mantle',
                                          'CN7 Epi Inner',
                                          'CN8 ASC Rich',
                                          'CN9 Connective') 

# Compute within-cluster cell count Z-scores across CNs
xp_l1l2clust_cn_mtx_scaled_t <- t(xp_l1l2clust_cn_mtx_scaled)

# Assemble heatmap
pdf('/filepath/fig4/fig4_supp/gc_closed_cn_l1_clust_distr_heatmap.pdf', 
    height = 5, width = 10,
    family = "sans")
cn_orig_l1_clust_distr_heatmap <- pheatmap(
  xp_l1l2clust_cn_mtx_scaled_t,
  color         = cn_heatmap_cols,
  breaks        = cn_heatmap_col_breaks,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  fontsize      = 13,
  fontsize_row  = 13,
  fontsize_col  = 13,
  angle_col     = "45",
  border_color  = "grey80",
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend        = TRUE
)
cn_orig_l1_clust_distr_heatmap
dev.off()

# Fig. S16I - Heatmap for L2 nnCD4 T cell subcluster enrichment across CNs (initial CN analysis version) ----

# Create matrix of L2 cluster cell counts per CN
xp_l1l2clust_cn_df <- xp_nncd4t_obj@meta.data %>%
  dplyr::select(cluster = l2_nncd4t_annot,
                neighborhood = cn_annot)
xp_l1l2clust_cn_table <- table(xp_l1l2clust_cn_df$cluster, xp_l1l2clust_cn_df$neighborhood)
xp_l1l2clust_cn_mtx <- as.matrix(xp_l1l2clust_cn_table)

# Compute within-cluster cell count Z-scores across CNs
xp_l1l2clust_cn_mtx_scaled <- t(scale(t(xp_l1l2clust_cn_mtx)))

# Build color palette for Z-scores centered around 0
cn_heatmap_col_breaks <- seq(-max(abs(xp_l1l2clust_cn_mtx_scaled)), max(abs(xp_l1l2clust_cn_mtx_scaled)), length.out = 100)
cn_heatmap_cols  <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(length(cn_heatmap_col_breaks) - 1))

# Specify CN names for heatmap
colnames(xp_l1l2clust_cn_mtx_scaled) <- c('CN0 Fol Border',
                                          'CN1 TCZ Inner',
                                          'CN2 GC LZ',
                                          'CN3 TCZ Outer',
                                          'CN4 Epi Outer',
                                          'CN5 GZ DZ',
                                          'CN6 Mantle',
                                          'CN7 Epi Inner',
                                          'CN8 ASC Rich',
                                          'CN9 Connective') 

# Transpose matrix to plot clusters as columns and CNs as rows on heatmap
xp_l1l2clust_cn_mtx_scaled_t <- t(xp_l1l2clust_cn_mtx_scaled)

# Assemble heatmap
pdf('/filepath/fig4/fig4_supp/cn_orig_l2nncd4t_clust_distr_heatmap.pdf', 
    height = 5, width = 10,
    family = "sans")
cn_orig_l1_clust_distr_heatmap <- pheatmap(
  xp_l1l2clust_cn_mtx_scaled_t,
  color         = cn_heatmap_cols,
  breaks        = cn_heatmap_col_breaks,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  fontsize      = 13,
  fontsize_row  = 13,
  fontsize_col  = 13,
  angle_col     = "45",
  border_color  = "grey80",
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend        = TRUE
)
cn_orig_l1_clust_distr_heatmap
dev.off()

# Fig. S16J - Heatmap for L2 nnCD4 T cell subcluster enrichment across CNs (merged GC CN analysis version) ----

# Create matrix of L2 cluster cell counts per CN
xp_l1l2clust_cn_df <- xp_nncd4t_obj@meta.data %>%
  dplyr::select(cluster = l2_nncd4t_annot,
                neighborhood = cn_gc_closed_annot)
xp_l1l2clust_cn_table <- table(xp_l1l2clust_cn_df$cluster, xp_l1l2clust_cn_df$neighborhood)
xp_l1l2clust_cn_mtx <- as.matrix(xp_l1l2clust_cn_table)

# Compute within-cluster cell count Z-scores across CNs
xp_l1l2clust_cn_mtx_scaled <- t(scale(t(xp_l1l2clust_cn_mtx)))

# Build color palette for Z-scores centered around 0
cn_heatmap_col_breaks <- seq(-max(abs(xp_l1l2clust_cn_mtx_scaled)), max(abs(xp_l1l2clust_cn_mtx_scaled)), length.out = 100)
cn_heatmap_cols  <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(length(cn_heatmap_col_breaks) - 1))

# Specify CN names for heatmap
colnames(xp_l1l2clust_cn_mtx_scaled) <- c('CN0 Fol Border',
                                          'CN1 TCZ Inner',
                                          'CN2 GC',
                                          'CN3 TCZ Outer',
                                          'CN4 Epi Outer',
                                          'CN6 Mantle',
                                          'CN7 Epi Inner',
                                          'CN8 ASC Rich',
                                          'CN9 Connective') 

# Transpose matrix to plot clusters as columns and CNs as rows on heatmap
xp_l1l2clust_cn_mtx_scaled_t <- t(xp_l1l2clust_cn_mtx_scaled)

# Assemble heatmap
pdf('/filepath/fig4/fig4_supp/gc_closed_cn_l2nncd4t_clust_distr_heatmap.pdf', 
    height = 5, width = 10,
    family = "sans")
cn_orig_l1_clust_distr_heatmap <- pheatmap(
  xp_l1l2clust_cn_mtx_scaled_t,
  color         = cn_heatmap_cols,
  breaks        = cn_heatmap_col_breaks,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  fontsize      = 13,
  fontsize_row  = 13,
  fontsize_col  = 13,
  angle_col     = "45",
  border_color  = "grey80",
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend        = TRUE
)
cn_orig_l1_clust_distr_heatmap
dev.off()