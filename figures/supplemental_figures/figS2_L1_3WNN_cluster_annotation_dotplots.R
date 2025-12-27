# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S2
# Fig. S2 – Trimodal annotation of Level 1 mononuclear cell states in tonsil and peripheral blood TEAseq dataset

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

# Import Seurat object from TEAseq Data Preprocessing Step 13 (trimodal dimensionality reduction and L1 3WNN clustering of all tonsil and peripheral blood mononuclear cells, including Harmony integration across donors ----
l1_teaseq_obj <- readRDS('/filepath/step13_bulk_harmony/l1_teaseq_obj.rds')

# Set cluster order for annotation dotplots
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
wnn_dotplot_order <- c('CD4 Tn','CD4 Tn Ribo','CD4 Tm','Treg','Tfh-like','CD8 Tn','CD8 Tm/Inn','CD8 UTC','NK','NBC','MBC','GCB','ASC','Myl','Myl CD16')
Idents(l1_teaseq_obj) <- factor(
  Idents(l1_teaseq_obj),
  levels = wnn_dotplot_order
)

# Set color palette
puor_colors <- brewer.pal(11, "PuOr")

# A) L1 3WNN Cluster Marker DotPlot - ATAC  ----
DefaultAssay(l1_teaseq_obj) <- 'ACT' # using Signac ATAC GeneActivity assay
wnn_clust_dotplot_features_atac <- c('TCF7','CD3E','CD4','LEF1','IL7R','ICOS','MAF','RORC','IL2','CCR8','FOXP3','IL2RA','CXCR5','PDCD1','BCL6','TOX2','IL21','GNG4','CD8A','CD8B','PRF1','GZMB','ZNF683','XCL1','NCR2','KLRD1','NCR1','PAX5','CD19','MS4A1','CD38','CD80','CD86','PRDM1','JCHAIN','TNFRSF17','CSF2RA','CD300LB','CD14','FCGR3A')
l1_wnn_clust_marker_dotplot_atac <- DotPlot(l1_teaseq_obj, features = wnn_clust_dotplot_features_atac, 
                                            cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_clust_dotplot_features_atac) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l1_teaseq_obj$l1_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10], # Blue
    mid = "white", # White at zero
    high = puor_colors[2], # Orange
    midpoint = 0 # White at zero
  )
pdf("/filepath/fig1/fig1_supp/l1_wnn_clust_marker_dotplot_atac.pdf", width = 14, height = 5)
l1_wnn_clust_marker_dotplot_atac
dev.off()

# B) L1 3WNN Cluster Marker DotPlot - RNA ----
DefaultAssay(l1_teaseq_obj) <- 'RNA'
wnn_clust_dotplot_features_rna <- c('CD3E','CD4','LEF1','RPS10','RPL31','IL7R','ICOS','GATA3','RORA','FOXP3','IL2RA','IKZF2','CXCR5','PDCD1','BCL6','TOX2','IL21','GNG4','CD8A','CD8B','KLRK1','PRF1','GZMB','TRDC','ZNF683','KLRC2','NCR1','KLRF1','PAX5','MS4A1','CD80','CD86','CD38','RGS13','JCHAIN','MZB1','CSF2RA','ITGAX','CD14','FCGR3A')
l1_wnn_clust_marker_dotplot_rna <- DotPlot(l1_teaseq_obj, features = wnn_clust_dotplot_features_rna, 
                                           cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_clust_dotplot_features_rna) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l1_teaseq_obj$l1_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10], # Blue
    mid = "white", # White at zero
    high = puor_colors[2], # Orange
    midpoint = 0 # White at zero
  )
pdf("/filepath/fig1/fig1_supp/l1_wnn_clust_marker_dotplot_rna.pdf", width = 14, height = 5)
l1_wnn_clust_marker_dotplot_rna
dev.off()

# C) L1 3WNN Cluster Marker DotPlot - ADT ----
DefaultAssay(l1_teaseq_obj) <- 'ADT'
wnn_clust_dotplot_features_adt <- c('adt-CD2','adt-CD3','adt-CD4','adt-CD45RA','adt-CCR7','adt-CD62L','adt-CD127','adt-CD49f','adt-CD278','adt-CD161','adt-CD95','adt-CD194','adt-CD25','adt-CD185','adt-CD279','adt-CD69','adt-TIGIT','adt-CD200','adt-CD272','adt-CD57','adt-CD8','adt-CD244','adt-GPR56','adt-TCR.Vd2','adt-TCR.Va7.2','adt-CD94','adt-CD314','adt-CD56','adt-CD335','adt-CD19','adt-CD21','adt-IgM','adt-IgD','adt-CD38','adt-CD71','adt-CLEC12A','adt-CD172a','adt-CD1c','adt-CD14','adt-CD16')
wnn_clust_dotplot_features_adt_trim_prefix <- sub("^adt-", "", wnn_clust_dotplot_features_adt)
l1_wnn_clust_marker_dotplot_adt <- DotPlot(l1_teaseq_obj, features = wnn_clust_dotplot_features_adt, 
                                           cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_clust_dotplot_features_adt) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l1_teaseq_obj$l1_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
l1_wnn_clust_marker_dotplot_adt$data$features.plot <- factor(
  l1_wnn_clust_marker_dotplot_adt$data$features.plot,
  levels = wnn_clust_dotplot_features_adt,
  labels = wnn_clust_dotplot_features_adt_trim_prefix
)
pdf("/filepath/fig1/fig1_supp/l1_wnn_clust_marker_dotplot_adt.pdf", width = 14, height = 5)
l1_wnn_clust_marker_dotplot_adt
dev.off()