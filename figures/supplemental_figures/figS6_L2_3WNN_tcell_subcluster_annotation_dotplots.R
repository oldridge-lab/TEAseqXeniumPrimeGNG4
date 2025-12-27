# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S6
# Fig. S6 - Trimodal annotation of Level 2 Tfh versus nonTfh states within TEAseq dataset.

# Set up working environment ----

# Set working directory
setwd('/filepath/fig1/fig1_supp/')

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

# Set color palette
puor_colors <- brewer.pal(11, "PuOr")

# Import L2 3WNN TEAseq T Cell Seurat object (subclustering from Step 14 of T Cells extracted from L1 object) ----
l2_teaseq_tcell_obj <- readRDS('/filepath/step14_tcell_subcluster/l2_teaseq_tcell_obj.rds')

# Set dotplot order of L2 T Cell subclusters
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'
wnn_tcell_dotplot_order <- c("CD4 Tn 1", "CD4 Tn 2", "CD4 Tn 3", "CD4 Tn 4", "CD4 Tn 5", "CD4 Tn 6", "CD4 Tn 7","CD4 Tem", "Th17", "Treg RORgt", "cTreg", "CD4 Tcm/fh", "Tfh IL10", "Tfh GC", "CD4 Inn", "PLZF Inn", "gdT", "CD8 Tn", "CD8 Tm", "CD8 CTL", "CD8 UTC")
Idents(l2_teaseq_tcell_obj) <- factor(
  Idents(l2_teaseq_tcell_obj),
  levels = wnn_tcell_dotplot_order
)

# A) L2 T Cell 3WNN Cluster Marker DotPlot - ATAC ----
DefaultAssay(l2_teaseq_tcell_obj) <- 'ACT' # using Signac ATAC GeneActivity assay
wnn_tcell_clust_dotplot_features_atac <- c('CD4','BACH2','LEF1','SELL','CCR7','IL7R','CD27','BHLHE40','TNF','CCR6','RORC','CCL20','IL2RA','FOXP3','KLF2','ICOS','MAF','IL6R','FKBP5','ST8SIA1','PTPN13','CXCR5','BCL6','TOX2','PDCD1','CXCL13','GNG4','IL21','IL10','PRDM1','ZBTB16','CEBPD','EOMES','KLRD1','ZEB2','CD8A','CD8B','GZMB','PRF1','ZNF683')
wnn_tcell_clust_marker_dotplot_atac <- DotPlot(l2_teaseq_tcell_obj, features = wnn_tcell_clust_dotplot_features_atac,
                                               cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_tcell_clust_dotplot_features_atac) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l2_teaseq_tcell_obj$l2_tcell_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
wnn_tcell_clust_marker_dotplot_atac
ggsave("/filepath/fig1/fig1_supp/wnn_tcell_clust_marker_dotplot_atac.pdf", plot = wnn_tcell_clust_marker_dotplot_atac, width = 14, height = 5.75)

# B) L2 T Cell 3WNN Cluster Marker DotPlot - RNA ----
DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
wnn_tcell_clust_dotplot_features_rna <- c('CD4','BACH2','LEF1','SELL','CCR7','IL7R','ITGB1','BHLHE40','CCR6','RORC','IL2RA','FOXP3','KLF2','ICOS','IL6R','IL6ST','MAF','PHACTR2','TSHZ2','FKBP5','ST8SIA1','PTPN13','CXCR5','BCL6','TOX2','PDCD1','CXCL13','GNG4','IL21','IL10','PRDM1','ZBTB16','CEBPD','TRGC1','TRDC','CD8A','CD8B','GNLY','GZMH','ZNF683')
wnn_tcell_clust_marker_dotplot_rna <- DotPlot(l2_teaseq_tcell_obj, features = wnn_tcell_clust_dotplot_features_rna,
                                              cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_tcell_clust_dotplot_features_rna) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l2_teaseq_tcell_obj$tcell_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
wnn_tcell_clust_marker_dotplot_rna
ggsave("/filepath/fig1/fig1_supp/wnn_tcell_clust_marker_dotplot_rna.pdf", plot = wnn_tcell_clust_marker_dotplot_rna, width = 14, height = 5.75)

# C) L2 T Cell 3WNN Cluster Marker DotPlot - ADT ----
DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
wnn_tcell_clust_dotplot_features_adt <- c('adt-CD3','adt-CD4','adt-CD45RA','adt-CD31','adt-CCR7','adt-CD62L','adt-CD162','adt-CD27','adt-CD127','adt-CD49f','adt-CD161','adt-CD194','adt-CD25','adt-CD39','adt-CD45RO','adt-CD38','adt-CD278','adt-CD69','adt-CD11a','adt-TIGIT','adt-CD185','adt-CD279','adt-CD304','adt-CD71','adt-CD200','adt-CD272','adt-CD57','adt-CD35','adt-TCR.Va7.2','adt-TCR.Vd2','adt-GPR56','adt-CD94','adt-KLRG1','adt-CD8','adt-CD244','adt-CD314','adt-CD16','adt-CD56','adt-CD335','adt-CD183')
wnn_tcell_clust_dotplot_features_adt_trim_prefix <- sub("^adt-", "", wnn_tcell_clust_dotplot_features_adt)
wnn_tcell_clust_marker_dotplot_adt <- DotPlot(l2_teaseq_tcell_obj, features = wnn_tcell_clust_dotplot_features_adt,
                                              cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_tcell_clust_dotplot_features_adt) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l2_teaseq_tcell_obj$tcell_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
wnn_tcell_clust_marker_dotplot_adt$data$features.plot <- factor(
  wnn_tcell_clust_marker_dotplot_adt$data$features.plot,
  levels = wnn_tcell_clust_dotplot_features_adt,
  labels = wnn_tcell_clust_dotplot_features_adt_trim_prefix
)
wnn_tcell_clust_marker_dotplot_adt
ggsave("/filepath/fig1/fig1_supp/wnn_tcell_clust_marker_dotplot_adt.pdf", plot = wnn_tcell_clust_marker_dotplot_adt, width = 14, height = 5.75)