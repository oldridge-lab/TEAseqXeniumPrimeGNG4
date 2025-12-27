# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S9
# Fig. S9 - SCENIC and chromVAR analyses infer distinct patterns of gene regulatory factor activity between Tfh states.

# Set up working environment ----

# Set working directory
setwd('/filepath/fig2/fig2_supp/')

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Matrix) # 1.7-0
library(tidyselect) # 1.2.1
library(BiocManager) # 1.30.23
library(GenomeInfoDb) # 1.40.0
library(GenomicRanges) # 1.56.0
library(IRanges) # 2.38.0
library(Rsamtools) # 2.20.0
library(S4Vectors) # 0.42.0
library(BiocGenerics) # 0.50.0
library(biovizBase) # 1.52.0
library(GenomicFeatures) # 1.56.0
library(EnsDb.Hsapiens.v86) # 2.99.0
library(BSgenome.Hsapiens.UCSC.hg38) # 1.4.5
library(Signac) # 1.13.0
library(Seurat) # 5.1.0
library(SeuratObject) # 5.0.2
library(data.table) # 1.15.4
library(readr) # 2.1.5
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4
library(stringr) # 1.5.1
library(glmGamPoi) # 1.16.0
library(TFBSTools) # 1.42.0
library(JASPAR2020) # 0.99.1
library(motifmatchr) # 1.26.0
library(ggseqlogo) # 0.2
library(caret) # 6.0-94
library(presto) # 1.0.0
library(SoupX) # 1.6.2
library(scDblFinder) # 1.19.1
library(Nebulosa) # 1.14.0
library(clusterProfiler) # 4.12.6 
library(enrichplot) # 1.24.2
library(org.Hs.eg.db) # 3.19.1
library(harmony) # 1.2.0
library(SCpubr) # 2.0.2
library(paletteer) # 1.6.0
library(RColorBrewer) # 1.1-3
library(EnhancedVolcano) # 1.22.0
library(AUCell) # 1.26.0
library(msigdbr) # 7.5.1
library(RcisTarget) # 1.23.1
library(GENIE3) # 1.26.0
library(SCENIC) # 1.3.1
library(arrow) # 18.1.0.1
library(doMC) # 1.3.8
library(R2HTML) # 2.3.4
library(reshape2) # 1.4.4
library(ggtern) # 3.5.0
library(ggh4x) # 0.3.0
library(SeuratWrappers) # 0.3.5
library(monocle3) # 1.3.7
library(cisTopic) # 0.3.0
library(densityClust) # 0.3.3
library(Rtsne) # 0.17
library(scatterplot3d) # 1.2
library(fastcluster) # 1.2.6
library(rtracklayer) # 1.64.0
library(tidyr) # 1.3.1
library(tidyverse) # 2.0.0
library(pheatmap) # 1.0.12
library(ggraph) # 2.2.1
library(igraph) # 2.0.3
library(BPCells) # 0.3.0
library(future) # 1.33.2
library(FNN) # 1.1.4
library(tibble) # 3.2.1
library(purrr) # 1.0.2
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # 3.18.0
library(ComplexHeatmap) # 2.20.0
library(circlize) # 0.4.16
library(BiocParallel) # 1.38.0
library(ggnewscale) # 0.5.0
library(ggVennDiagram) # 1.5.2
library(biomaRt) # 2.60.0
library(writexl) # 1.5.4

# Import Seurat object from TEAseq Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors) ----
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Assign names to Tfh subcluster numbers
tfh_subcluster_names <- c(
  "0" = "Tfh-Circ",
  "1" = "Tcm",
  "2" = "Tfh-Int",
  "3" = "Tfh-NFATC1",
  "4" = "Tfh-IL10",
  "5" = "Tfh-Resting",
  "6" = "Tfh-CXCL13",
  "7" = "Tfh-BOB1",
  "8" = "Tfh-AP1"
)

# Specify object idents
l3_teaseq_tfh_obj$tfh_wnn_annot <- l3_teaseq_tfh_obj$wnn.tfh.subcluster.harmony
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
l3_teaseq_tfh_obj <- RenameIdents(l3_teaseq_tfh_obj, tfh_subcluster_names)
l3_teaseq_tfh_obj$tfh_wnn_annot <- Idents(l3_teaseq_tfh_obj)
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
table(Idents(l3_teaseq_tfh_obj))

# A) Heatmap of SCENIC regulon AUC scores scaled across L3 Tfh-like clusters ----
setwd('/filepath/l3_tfh_scenic_output_folder/')
tfh_scenicOptions <- readRDS("/filepath/l3_tfh_scenic_output_folder/scenicOptions4.rds") # complete run through binarization step
tfh_regulonAUC <- loadInt(tfh_scenicOptions, "aucell_regulonAUC") # load AUC regulon matrix
tfh_AUCmatrix <- getAUC(tfh_regulonAUC) # retrieve AUC values
tfh_cellInfo <- tfh_scenicOptions@inputDatasetInfo$cellInfo # retrieve cell metadata - note SCENIC was run with preliminary cluster annotations that we later updated in final Preprocessing Step 15 code

# Scaled regulon activity by cluster (all regulons) - duplicates filtered
tfh_regulonAUC_filt <- tfh_regulonAUC[onlyNonDuplicatedExtended(rownames(tfh_regulonAUC)),] # removed 64 duplicated rownames - function 'returns the regulon names filtering-out the "extended" regulons if there is a regulon based on high-confidence annotations'
tfh_clust_regulonActivity <- sapply(split(rownames(tfh_cellInfo), tfh_cellInfo$tfh_wnn_annot),
                                    function(cells) rowMeans(getAUC(tfh_regulonAUC_filt)[,cells]))
tfh_clust_regulonActivity_scaled <- t(scale(t(tfh_clust_regulonActivity), center = T, scale=T))

# Set colors palette
puor_colors <- brewer.pal(11, "PuOr")
tfh_scenic_heatmap_min_val <- min(tfh_clust_regulonActivity_scaled, na.rm = TRUE)
tfh_scenic_heatmap_max_val <- max(tfh_clust_regulonActivity_scaled, na.rm = TRUE)
tfh_scenic_heatmap_hi_col  <- brewer.pal(11, "PuOr")[2]
tfh_scenic_heatmap_lo_col <- brewer.pal(11, "PuOr")[10]
tfh_scenic_heatmap_cols  <- colorRamp2(c(tfh_scenic_heatmap_min_val, 0, tfh_scenic_heatmap_max_val),
                                       c(tfh_scenic_heatmap_lo_col, "white", tfh_scenic_heatmap_hi_col))

# Set cluster order for heatmap
wnn_tfh_heatmap_order <- c("Tcm", "Tfh-Circ", "Tfh-AP1", "Tfh-Resting", "Tfh-IL10", "Tfh-Int", "Tfh-NFATC1","Tfh-CXCL13", "Tfh-BOB1")
tfh_clust_regulonActivity_scaled <- tfh_clust_regulonActivity_scaled[, wnn_tfh_heatmap_order]

# Get activity of top regulons per cluster
topRegulators <- melt(tfh_clust_regulonActivity_scaled)
colnames(topRegulators) <- c("Regulon", "tfh_subset", "RelativeActivity")
topRegulators <- topRegulators %>% # filter for top five active regulons
  filter(RelativeActivity > 0.5) %>%
  group_by(tfh_subset) %>%
  slice_max(order_by = RelativeActivity, n = 9) %>%
  ungroup()
topRegulonNames <- unique(topRegulators$Regulon) # get union of regulons selected across clusters
regulonActivity_top <- tfh_clust_regulonActivity_scaled[rownames(tfh_clust_regulonActivity_scaled) %in% topRegulonNames, ] # subset original scaled matrix to include only top regulons

# Trim regulon names to TF root for visualization
rownames(regulonActivity_top) <- ifelse(
  grepl("_extended", rownames(regulonActivity_top)),
  sub("_extended.*", "", rownames(regulonActivity_top)),          
  sub(" .*", "", rownames(regulonActivity_top))                     
)

# Assemble SCENIC heatmap
pdf('/filepath/fig2/fig2_supp/tfh_scenic_chromvar_heatmaps/tfh_scenic_supplement_heatmap_test.pdf', width = 20, height = 5.75)
draw(Heatmap(
  t(regulonActivity_top),
  name = NULL,                          
  col = tfh_scenic_heatmap_cols,                         
  show_row_dend = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45, 
  column_names_gp = gpar(fontsize = 12, fontfamily = "sans"),
  row_names_side = 'left',
  row_names_gp = gpar(fontsize = 15, fontfamily = "sans"),
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    labels_gp = gpar(fontsize = 12, fontfamily = "sans"))
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# B) Heatmap of chromVAR TF motif accessibility deviation scores scaled across L3 Tfh-like clusters (ATAC-based inference) ----

# Select chromVAR TF of interest enriched in any cluster from FindAllMarkers analysis with pvaladj < 0.05
l3_chromvar_heatmap_tf_list <- c('HIF1A','CTCF','YY1','RFX5','GLIS3','E2F2','VEZF1', # Tcm
                                 'RFX2','RFX7','Klf1', # Tfh-Circ
                                 'ELK1','SP1','KLF2','KLF3','SP3','KLF14','FEV','ETV3', # Tfh-AP1
                                 'SP4','KLF6','ETV2','KLF13','GABPA','STAT1','Klf12','ERF', # Tfh-Resting
                                 'IRF9','IRF4','MAFK','NFE2L1','IRF8','STAT3', # Tfh-IL10
                                 'NR3C2','NR3C1','Ar','POU5F1B','POU3F4','ZSCAN4','REL','ZBTB32', # Tfh-Int
                                 'NFATC3','NFATC1','NFATC4','NFAT5','NFATC2', # Tfh-NFATC1
                                 'BATF3','BATF','BACH1','JUNB','JDP2','CEBPA','BACH2','FOSL1','FOSL2', # Tfh-CXCL13
                                 'POU5F1','POU2F3','Ascl2','MYOG','Tcf12','MAF','IKZF1','ASCL1','ZBTB18' # Tfh-BOB1
)
tfh_chromvar_mtx_full <- GetAssayData(l3_teaseq_tfh_obj, assay = 'chromvar', slot = "data")

# Map chromVAR motif names to TF names for heatmap visualization
pfm <- getMatrixSet( 
  x = JASPAR2020, # get TF motif PFM information
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
id_to_tf <- sapply(pfm, function(x) name(x))
names(id_to_tf) <- sapply(pfm, ID)
rn_tfh_chromvar <- rownames(tfh_chromvar_mtx_full)
rn_tfh_chromvar_tfname <- id_to_tf[rn_tfh_chromvar]
rownames(tfh_chromvar_mtx_full) <- rn_tfh_chromvar_tfname

# Prepare heatmap data matrix
chromvar_tf_list <- intersect(rownames(tfh_chromvar_mtx_full), l3_chromvar_heatmap_tf_list) # subset matrix to enriched TF of interest
tfh_chromvar_mtx <- as.matrix(tfh_chromvar_mtx_full[chromvar_tf_list, , drop = FALSE])  # features x cells
dim(tfh_chromvar_mtx)
tfh_idents <- Idents(l3_teaseq_tfh_obj)
tfh_ident_levels <- levels(tfh_idents)
tfh_chromvar_mtx_by_cluster <- sapply(tfh_ident_levels, function(k) {
  idx <- which(tfh_idents == k)
  rowMeans(tfh_chromvar_mtx[, idx, drop = FALSE], na.rm = TRUE)
})
tfh_chromvar_mtx_t <- t(tfh_chromvar_mtx_by_cluster) # transpose matrix to be clusters x TF names
colnames(tfh_chromvar_mtx_t) <- rownames(tfh_chromvar_mtx)
rownames(tfh_chromvar_mtx_t) <- tfh_ident_levels
tfh_chromvar_mtx_scaled <- as.matrix(t(scale(t(tfh_chromvar_mtx_t), center = TRUE, scale = TRUE)))

# Set colors
tfh_chromvar_heatmap_min_val <- min(tfh_chromvar_mtx_scaled, na.rm = TRUE)
tfh_chromvar_heatmap_max_val <- max(tfh_chromvar_mtx_scaled, na.rm = TRUE)
tfh_chromvar_heatmap_hi_col  <- brewer.pal(11, "PuOr")[2]
tfh_chromvar_heatmap_lo_col <- brewer.pal(11, "PuOr")[10]
tfh_chromvar_heatmap_cols  <- colorRamp2(c(tfh_chromvar_heatmap_min_val, 0, tfh_chromvar_heatmap_max_val),
                                         c(tfh_chromvar_heatmap_lo_col, "white", tfh_chromvar_heatmap_hi_col))

# Assemble chromVAR TF heatmap
pdf('/filepath/fig2/fig2_supp/tfh_scenic_chromvar_heatmaps/tfh_chromvar_supplement_heatmap.pdf', width = 16, height = 5.75)
draw(Heatmap(
  tfh_chromvar_mtx_scaled,
  name = NULL,                          
  col = tfh_chromvar_heatmap_cols,                         
  show_row_dend = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45, 
  column_names_gp = gpar(fontsize = 12, fontfamily = "sans"),
  row_names_side = 'left',
  row_names_gp = gpar(fontsize = 15, fontfamily = "sans"),
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    labels_gp = gpar(fontsize = 12, fontfamily = "sans"))
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()