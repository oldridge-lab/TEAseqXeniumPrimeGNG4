# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S8
# Fig. S8 - Trimodal annotation of Level 3 Tfh versus Tcm states within TEAseq dataset.

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

# Import Seurat object from TEAseq Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors)
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

# Set color palette
puor_colors <- brewer.pal(11, "PuOr")

# A) L3 Tfh 3WNN Cluster Marker DotPlot - ATAC ----
DefaultAssay(l3_teaseq_tfh_obj) <- 'ACT' # (using Signac ATAC GeneActivity assay)

# Select features of interest
wnn_tfh_clust_dotplot_features_atac <- c('CD55','KLF2','ANXA1','S1PR1','IL2RA','CCR7','SELPLG','SELL',
                                         'ADAM8','SH2D3A','TGFB1','SKI',
                                         'CD96','DOCK2','MAFK','KCNK4',
                                         'IRF7','IL7R','IKZF1',
                                         'FOXB1','PRDM1','IL10','CTLA4',
                                         'GNG4','PTPN14','SMCO4','XXYLT1','MSI2','PLCL2',
                                         'NFATC1','TOX','PTPN11','BTLA','CD200',
                                         'CXCL13','IL12RB2','ASCL2','BCAT1','SH2D1A','CD38',
                                         'POU2AF1','PDCD1','CD40LG','TOX2','ICOS','B3GAT1','IL21','CXCR5','TIGIT','ST8SIA1'
)
wnn_tfh_clust_marker_dotplot_atac <- DotPlot(l3_teaseq_tfh_obj, features = wnn_tfh_clust_dotplot_features_atac,
                                             cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_tfh_clust_dotplot_features_atac) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l3_teaseq_tfh_obj$tfh_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
wnn_tfh_clust_marker_dotplot_atac
pdf("/filepath/fig2/fig2_supp/wnn_tfh_clust_marker_dotplot_atac.pdf", width = 17, height = 3.75)
wnn_tfh_clust_marker_dotplot_atac
dev.off()

# B) L3 Tfh 3WNN Cluster Marker DotPlot - RNA ----
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
wnn_tfh_clust_dotplot_features_rna <- c('KLF2','SELL','CCR7','ANXA1','CD55',
                                        'ANK3','YPEL5','CAMK4','IL7R','FOXP1','SARAF',
                                        'JUNB','JUN','FOSB','NR4A1','NFKBIA',
                                        'PDE7A','FYB1','THEMIS','IKZF1','CD40LG',
                                        'PRDM1','IL10','IL2RA','LAG3','MAF','CD38','IKZF3',
                                        'TOX2','TOX','DAB1','KSR2','PTPN14',
                                        'NFATC1','CD200','XXYLT1','BTLA','PDCD1',
                                        'CXCL13','IL21','CTLA4','SH2D1A','CXCR5',
                                        'POU2AF1','B3GAT1','DRAIC','GNG4','PLCG2','IKZF2','TIGIT'
)
wnn_tfh_clust_marker_dotplot_rna <- DotPlot(l3_teaseq_tfh_obj, features = wnn_tfh_clust_dotplot_features_rna,
                                            cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_tfh_clust_dotplot_features_rna) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l3_teaseq_tfh_obj$tfh_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
wnn_tfh_clust_marker_dotplot_rna
pdf("/filepath/fig2/fig2_supp/wnn_tfh_clust_marker_dotplot_rna.pdf", width = 17, height = 3.75)
wnn_tfh_clust_marker_dotplot_rna
dev.off()

# C) L3 Tfh 3WNN Cluster Marker DotPlot - ADT ----
DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
wnn_tfh_clust_dotplot_features_adt <- c('adt-CD2','adt-CD3','adt-CD4','adt-CD27','adt-CD49d','adt-CD29','adt-CD44','adt-CD162','adt-CD62L','adt-CCR7','adt-CD26','adt-CD352',
                                        'adt-CD48','adt-CD127','adt-CD5','adt-CD55','adt-CD196',
                                        'adt-CD49f','adt-CD99',
                                        'adt-CD7','adt-CD49a','adt-CD47','adt-CD183',
                                        'adt-CD161','adt-CD25','adt-CD81','adt-CD195','adt-CD45RO','adt-CD39','adt-CD151','adt-CD194','adt-CD38',
                                        'adt-CD304','adt-CD305','adt-CD278','adt-CD279','adt-CD11a','adt-CD69','adt-HLA.DR','adt-CD85j','adt-CD18',
                                        'adt-CD71','adt-CD57','adt-CD154','adt-CD134','adt-CD185','adt-CD192',
                                        'adt-TIGIT','adt-CD58','adt-CD272'
)
wnn_tfh_clust_dotplot_features_adt_trim_prefix <- sub("^adt-", "", wnn_tfh_clust_dotplot_features_adt)
wnn_tfh_clust_marker_dotplot_adt <- DotPlot(l3_teaseq_tfh_obj, features = wnn_tfh_clust_dotplot_features_adt,
                                            cluster.idents = FALSE, cols = 'PuOr') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain', size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(wnn_tfh_clust_dotplot_features_adt) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l3_teaseq_tfh_obj$tfh_wnn_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
wnn_tfh_clust_marker_dotplot_adt$data$features.plot <- factor(
  wnn_tfh_clust_marker_dotplot_adt$data$features.plot,
  levels = wnn_tfh_clust_dotplot_features_adt,
  labels = wnn_tfh_clust_dotplot_features_adt_trim_prefix
)
wnn_tfh_clust_marker_dotplot_adt
pdf("/filepath/fig2/fig2_supp/wnn_tfh_clust_marker_dotplot_adt.pdf", width = 17, height = 3.75)
wnn_tfh_clust_marker_dotplot_adt
dev.off()