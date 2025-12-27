# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S5
# Fig. S5 - Separate unimodal analyses of trimodal TEAseq data identify shared and distinct immune cell states. 

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

# Specify colors for tissues
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# Import Seurat object from Step 13 (RNA ATAC ADT 3WNN dimensionality reduction and L1 clustering of all PBMC & Tonsil mononuclear cells including Harmony integration across donors ----
l1_teaseq_obj <- readRDS('/filepath/l1_teaseq_obj.rds')

# Unimodal ATAC - L1 DR & Clustering ----

# Unimodal L1 ATAC cluster resolution for 15x clusters
DefaultAssay(l1_teaseq_obj) <- "ATAC"
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, resolution = 0.45, cluster.name = "atac.bulk.cluster.harmony", assay = 'ATAC', graph.name = 'ATAC_harmony_snn') # 0.45 Res for 15 Clusters
Idents(l1_teaseq_obj) <- "atac.bulk.cluster.harmony"

# Level 1 ATAC Cluster Annotation (15x clusters)
l1_atac_clust_names <- c(
  "0" = "CD4 Tn 1",
  "1" = "CD4 Tn 2",
  "2" = "Tfh-like",
  "3" = "NBC 1",
  "4" = "CD4 Tm 1",
  "5" = "CD8 Tm/Inn",
  "6" = "CD4 Tn 3",
  "7" = "MBC",
  "8" = "CD8 Tn",
  "9" = "NK",
  "10" = "GCB",
  "11" = "Myl",
  "12" = "NBC 2",
  "13" = "CD4 Tm 2",
  "14" = "ASC"
)

# Rename unimodal ATAC L1 15x cluster numbers to assigned names
l1_teaseq_obj$l1_atac_annot <- l1_teaseq_obj$atac.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_atac_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_atac_clust_names)
l1_teaseq_obj$l1_atac_annot <- Idents(l1_teaseq_obj)

# Specify L1 unimodal ATAC cluster colors
l1_atac_clust_cols <- c(
  "CD4 Tn 1" = "#8ab1cc",       
  "CD4 Tn 2" = "#535c81",
  "Tfh-like" = "#e49a78",         
  "NBC 1" = "#657ab0",        
  "CD4 Tm 1" = "#b37cbf",      
  "CD8 Tm/Inn" = "#add5df",
  "CD4 Tn 3" = "#c9744d",
  "MBC" = "#d4b5e4",         
  "CD8 Tn" = "#50664f",      
  "NK" = "#728782",
  "GCB" = "#d66d97",         
  "Myl" = "#b6b7ba",         
  "NBC 2" = "#f4a6c7",
  "CD4 Tm 2" = "#72d1b4",        
  "ASC" = "#dfcc78"          
)

# Unimodal ATAC L1 UMAP
l1_atac_umap_atac_clust <- DimPlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = TRUE, label.box = TRUE, label.size = 6, pt.size = 0.75, repel = TRUE, group.by = 'l1_atac_annot', cols = l1_atac_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_atac_umap_atac_clust
ggsave("/filepath/fig1/fig1_supp/l1_atac_unimodal/l1_atac_umap_atac_clust.pdf", plot = l1_atac_umap_atac_clust, width = 6, height = 6)

# Unimodal ATAC L1 Cluster Annotation DotPlot
puor_colors <- brewer.pal(11, "PuOr")
Idents(l1_teaseq_obj) <- 'l1_atac_annot'
atac_dotplot_order <- c('CD4 Tn 1','CD4 Tn 2','CD4 Tn 3','CD4 Tm 1','CD4 Tm 2','Tfh-like','CD8 Tn','CD8 Tm/Inn','NK','NBC 1','NBC 2','MBC','GCB','ASC','Myl')
Idents(l1_teaseq_obj) <- factor(
  Idents(l1_teaseq_obj),
  levels = atac_dotplot_order
)
DefaultAssay(l1_teaseq_obj) <- 'ACT' # Using Signac GeneActivity assay
atac_dotplot_lin_features <- c('TCF7','CD3E','CD4','LEF1','IL7R','FOXP3','IL2RA','RORC','IL2','ICOS','MAF','CTLA4','CD200','CXCR5','PDCD1','BCL6','TOX2','IL21','GNG4','CD8A','CD8B','PRF1','KLRK1','GZMB','ZNF683','XCL1','KLRD1','NCR1','PAX5','MS4A1','CD80','CD86','CD38','JCHAIN','TNFRSF17','PRDM1','CSF2RA','CD300LB','CD14','FCGR3A')
l1_atac_clust_annot_dotplot <- DotPlot(l1_teaseq_obj, features = atac_dotplot_lin_features, 
                                       cluster.idents = FALSE, cols = 'PuOr', dot.scale = 5.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(atac_dotplot_lin_features) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l1_teaseq_obj$atac.bulk.cluster.harmony))    - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10], # Blue
    mid = "white", # White at zero
    high = puor_colors[2], # Orange
    midpoint = 0 # White at zero
  )
l1_atac_clust_annot_dotplot
ggsave("/filepath/fig1/fig1_supp/l1_atac_unimodal/l1_atac_clust_annot_dotplot.pdf", plot = l1_atac_clust_annot_dotplot, width = 14, height = 5)

# Flipping dotplot for visualization
l1_atac_clust_annot_dotplot_flip <- l1_atac_clust_annot_dotplot + 
  coord_flip() +
  theme(
    axis.text.x = element_text(face = 'plain'),
    axis.text.y = element_text(face = 'italic'),
    legend.position = "bottom",       
    legend.direction = "horizontal",  
    legend.box = "horizontal",         
    legend.justification = "center",
    legend.margin = margin(t = 5, r = 5, b = 5, l = -20),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
  ) +
  guides(
    size = guide_legend(order = 1,  title.position = "top", title = 'Percent Expressed', title.hjust = 0.5),
    color = guide_colourbar(order = 2, title.position = "top", title = 'Average Expression', title.hjust = 0.5))
ggsave("/filepath/fig1/fig1_supp/l1_atac_unimodal/l1_atac_clust_annot_dotplot_flip.pdf", plot = l1_atac_clust_annot_dotplot_flip, width = 5, height = 10)

# Unimodal RNA - L1 DR & Clustering ----

# Unimodal L1 RNA cluster resolution for 15x clusters
DefaultAssay(l1_teaseq_obj) <- "RNA"
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, assay = 'RNA', resolution = 0.45, cluster.name = "rna.bulk.cluster.harmony") # 0.45 Res for 15 Clusters
Idents(l1_teaseq_obj) <- "rna.bulk.cluster.harmony"

# L1 RNA Cluster Annotation (15x clusters)
l1_rna_clust_names <- c(
  "0" = "CD4 Tn",
  "1" = "CD4 Tm",
  "2" = "NBC",
  "3" = "CD8 Tm/Inn 1",
  "4" = "Tfh-like",
  "5" = "MBC",
  "6" = "CD8 Tn",
  "7" = "Treg",
  "8" = "GCB",
  "9" = "NK",
  "10" = "CD8 Tm/Inn 2",
  "11" = "Myl",
  "12" = "Myl CD16",
  "13" = "CD4 Tn Ribo",
  "14" = "ASC"
)

# Rename L1 RNA 15x cluster numbers to assigned names
l1_teaseq_obj$l1_rna_annot <- l1_teaseq_obj$rna.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_rna_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_rna_clust_names)
l1_teaseq_obj$l1_rna_annot <- Idents(l1_teaseq_obj)

# Specify colors for unimodal L1 RNA clusters
l1_rna_clust_cols <- c(
  "CD4 Tn" = "#8ab1cc",       
  "CD4 Tm" = "#94a890",     
  "NBC" = "#657ab0",        
  "CD8 Tm/Inn 1" = "#b37cbf",
  "Tfh-like" = "#e49a78",      
  "MBC" = "#d4b5e4",         
  "CD8 Tn" = "#c49980",      
  "Treg" = "#c9a3c1",        
  "GCB" = "#d66d97",         
  "NK" = "#728782",          
  "CD8 Tm/Inn 2" = "#add5df",
  "Myl" = "#b6b7ba",         
  "Myl CD16" = "#e1d4c7",    
  "CD4 Tn Ribo" = "#e0c1b4",        
  "ASC" = "#dfcc78"          
)

# Unimodal L1 RNA UMAP
rotate_2D <- function(X, theta) { # function to rotate UMAP coordinate matrix, does not change underlying local or global cell-cell UMAP distances, purely for visualization purposes
  mat_2D_rot <- matrix(
    c(cos(theta), sin(theta), -sin(theta), cos(theta)),
    2,
    2
  )
  Xp <- t(mat_2D_rot %*% t(X))
  return(Xp)
}
l1_teaseq_obj@reductions$umap.rna.harmony_rot <- l1_teaseq_obj@reductions$umap.rna.harmony
l1_teaseq_obj@reductions$umap.rna.harmony_rot@cell.embeddings <- rotate_2D(l1_teaseq_obj@reductions$umap.rna.harmony_rot@cell.embeddings, pi)
colnames(l1_teaseq_obj@reductions$umap.rna.harmony_rot@cell.embeddings) <- colnames(l1_teaseq_obj@reductions$umap.rna.harmony@cell.embeddings)
l1_rna_umap_rna_clust_rot <- DimPlot(l1_teaseq_obj, reduction = "umap.rna.harmony_rot", label = TRUE, label.box = TRUE, label.size = 6, pt.size = 0.75, repel = TRUE, group.by = 'l1_rna_annot', cols = l1_rna_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_rna_umap_rna_clust_rot
ggsave("/filepath/fig1/fig1_supp/l1_rna_unimodal/l1_rna_umap_rna_clust_rot.pdf", plot = l1_rna_umap_rna_clust_rot, width = 6, height = 6)

# Unimodal L1 RNA Cluster Annotation DotPlot
puor_colors <- brewer.pal(11, "PuOr")
DefaultAssay(l1_teaseq_obj) <- 'RNA'
rna_dotplot_order <- c('CD4 Tn','CD4 Tn Ribo','CD4 Tm','Treg','Tfh-like','CD8 Tn','CD8 Tm/Inn 1','CD8 Tm/Inn 2','NK','NBC','MBC','GCB','ASC','Myl','Myl CD16')
Idents(l1_teaseq_obj) <- factor(
  Idents(l1_teaseq_obj),
  levels = rna_dotplot_order
)
rna_dotplot_lin_features <- c('TCF7','CD3E','CD4','CCR7','IL7R','RPS27','RPL14','RORA','GATA3','ICOS','FOXP3','IL2RA','CXCR5','PDCD1','BCL6','TOX2','IL21','GNG4','CD8A','CD8B','KLRK1','ZNF683','TRDC','PRF1','RUNX3','GZMB','KLRC2','NCR1','KLRF1','MS4A1','IGHM','CD80','RGS13','XBP1','JCHAIN','MZB1','CSF2RA','ITGAX','CD14','FCGR3A')
l1_rna_clust_annot_dotplot <- DotPlot(l1_teaseq_obj, features = rna_dotplot_lin_features, 
                                      cluster.idents = FALSE, cols = 'PuOr', dot.scale = 5.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(rna_dotplot_lin_features) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l1_teaseq_obj$rna.bulk.cluster.harmony))    - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10], # Blue
    mid = "white", # White at zero
    high = puor_colors[2], # Orange
    midpoint = 0 # White at zero
  )
l1_rna_clust_annot_dotplot
ggsave("/filepath/fig1/fig1_supp/l1_rna_unimodal/l1_rna_clust_annot_dotplot.pdf", plot = l1_rna_clust_annot_dotplot, width = 12, height = 6)

# Flipping DotPlot for visualization
l1_rna_clust_annot_dotplot_flip <- l1_rna_clust_annot_dotplot + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain'),
    axis.text.y = element_text(face = 'italic'),
    legend.position = "bottom",       
    legend.direction = "horizontal",  
    legend.box = "horizontal",         
    legend.justification = "center",    
    legend.margin = margin(t = 5, r = 5, b = 5, l = -20),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
  ) +
  guides(
    size = guide_legend(order = 1,  title.position = "top", title = 'Percent Expressed', title.hjust = 0.5),
    color = guide_colourbar(order = 2, title.position = "top", title = 'Average Expression', title.hjust = 0.5))
ggsave("/filepath/fig1/fig1_supp/l1_rna_unimodal/l1_rna_clust_annot_dotplot_flip.pdf", plot = l1_rna_clust_annot_dotplot_flip, width = 5, height = 10)

# Unimodal ADT - L1 DR & Clustering ----

# Unimodal L1 ADT cluster resolution for 15x clusters
DefaultAssay(l1_teaseq_obj) <- "ADT"
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, resolution = 0.84, cluster.name = "adt.bulk.cluster.harmony") # 0.84 Res for 15 Clusters
Idents(l1_teaseq_obj) <- "adt.bulk.cluster.harmony"

# Unimodal L1 ADT cluster annotation (15x clusters)
l1_adt_clust_names <- c(
  "0" = "CD4 Tn 1",
  "1" = "NBC",
  "2" = "CD4 Tm",
  "3" = "Tfh-like",
  "4" = "CD8 Tm/Inn 1",
  "5" = "Ag Exp B",
  "6" = "CD8 Tn",
  "7" = "NK",
  "8" = "CD8 Tm/Inn 2",
  "9" = "CD4 Tn 2",
  "10" = "Mono CD14",
  "11" = "Mono CD16",
  "12" = "CD8 Tm/Inn 3",
  "13" = "DC",
  "14" = "MAIT"
)

# Rename Level 1 coarse 15 ADT cluster numbers to assigned names
l1_teaseq_obj$l1_adt_annot <- l1_teaseq_obj$adt.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_adt_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_adt_clust_names)
l1_teaseq_obj$l1_adt_annot <- Idents(l1_teaseq_obj)

# Colors for unimodal L1 ADT clusters
l1_adt_clust_cols <- c(
  "CD4 Tn 1" = "#8ab1cc",       
  "NBC" = "#90a6bf",        
  "CD4 Tm" = "#94a890",      
  "Tfh-like" = "#e49a78",         
  "CD8 Tm/Inn 1" = "#add5df",
  "Ag Exp B" = "#d66d97",     
  "CD8 Tn" = "#c49980",   
  "NK" = "#728782",          
  "CD8 Tm/Inn 2" = "#b6b7ba",
  "CD4 Tn 2" = "#f4a6c7",
  "Mono CD14" = "#535c81", 
  "Mono CD16" = "#e1d4c7",    
  "CD8 Tm/Inn 3" = "#b37cbf",        
  "DC" = "#c9744d",     
  "MAIT" = "#72d1b4"     
)

# Unimodal L1 ADT UMAP
l1_adt_umap_adt_clust <- DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = TRUE, label.box = TRUE, label.size = 6, pt.size = 0.75, repel = TRUE, group.by = 'l1_adt_annot', cols = l1_adt_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_adt_umap_adt_clust
ggsave("/filepath/fig1/fig1_supp/l1_adt_unimodal/l1_adt_umap_adt_clust.pdf", plot = l1_adt_umap_adt_clust, width = 6, height = 6)

# Unimodal L1 ADT Cluster Annotation DotPlot
DefaultAssay(l1_teaseq_obj) <- 'ADT'
Idents(l1_teaseq_obj) <- 'l1_adt_annot'
adt_dotplot_order <- c('CD4 Tn 1','CD4 Tn 2','CD4 Tm','MAIT','Treg','Tfh-like','CD8 Tn','CD8 Tm/Inn 1','CD8 Tm/Inn 2','CD8 Tm/Inn 3','NK','NBC','Ag Exp B','Mono CD14','Mono CD16','DC')
Idents(l1_teaseq_obj) <- factor(
  Idents(l1_teaseq_obj),
  levels = adt_dotplot_order
)
adt_dotplot_lin_features <- c('adt-CD2','adt-CD3','adt-CD4','adt-CD27','adt-CD45RA','adt-CD62L','adt-CD127','adt-CCR7','adt-CD49f','adt-CD103','adt-CD95','adt-CD194','adt-CD25','adt-TCR.Va7.2','adt-CD278','adt-CD185','adt-CD279','adt-CD69','adt-TIGIT','adt-CD200','adt-CD8','adt-CD314','adt-CD94','adt-GPR56','adt-TCR.Vd2','adt-CD158e1','adt-CD56','adt-CD244','adt-CD335','adt-CD19','adt-CD21','adt-IgD','adt-CD38','adt-CD71','adt-CLEC12A','adt-CD172a','adt-CD14','adt-CD16','adt-HLA.DR.DP.DQ','adt-CD1c')
adt_dotplot_lin_features_trim_prefix <- sub("^adt-", "", adt_dotplot_lin_features)
l1_adt_clust_annot_dotplot <- DotPlot(l1_teaseq_obj, features = adt_dotplot_lin_features, 
                                      cluster.idents = FALSE, cols = 'PuOr', dot.scale = 5.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = seq(1.5, length(adt_dotplot_lin_features) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l1_teaseq_obj$adt.bulk.cluster.harmony))    - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  scale_color_gradient2(
    low = puor_colors[10], # Blue
    mid = "white", # White at zero
    high = puor_colors[2], # Orange
    midpoint = 0 # White at zero
  )
l1_adt_clust_annot_dotplot$data$features.plot <- factor(
  l1_adt_clust_annot_dotplot$data$features.plot,
  levels = adt_dotplot_lin_features,
  labels = adt_dotplot_lin_features_trim_prefix
)
ggsave("/filepath/fig1/fig1_supp/l1_adt_unimodal/l1_adt_clust_annot_dotplot.pdf", plot = l1_adt_clust_annot_dotplot, width = 12, height = 6)

# Rotate DotPlot for visualization
l1_adt_clust_annot_dotplot_flip <- l1_adt_clust_annot_dotplot + 
  coord_flip() +
  theme(
    legend.position = "bottom",       
    legend.direction = "horizontal",  
    legend.box = "horizontal",         
    legend.justification = "center",    
    legend.margin = margin(t = 5, r = 5, b = 5, l = -20),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
  ) +
  guides(
    size = guide_legend(order = 1,  title.position = "top", title = 'Percent Expressed', title.hjust = 0.5),
    color = guide_colourbar(order = 2, title.position = "top", title = 'Average Expression', title.hjust = 0.5))
ggsave("/filepath/fig1/fig1_supp/l1_adt_unimodal/l1_adt_clust_annot_dotplot_flip.pdf", plot = l1_adt_clust_annot_dotplot_flip, width = 5, height = 10)

# Determining features used in DotPlots above to annotate L1 unimodal clusters ----

# TEAseq Unimodal ATAC (Peaks) L1 Cluster Markers
DefaultAssay(l1_teaseq_obj) <- 'ATAC'
Idents(l1_teaseq_obj) <- 'atac.bulk.cluster.harmony'
l1_teaseq_unimodal_atac_clust_dap <- FindAllMarkers(l1_teaseq_obj, assay = 'ATAC') # default parameters
l1_teaseq_unimodal_atac_clust_dap$l1_unimodal_atac_clust_annot <- l1_atac_clust_names[as.character(l1_teaseq_unimodal_atac_clust_dap$cluster)]
dap_names <- l1_teaseq_unimodal_atac_clust_dap$gene
dap_gr <- StringToGRanges(dap_names)
dap_closest_gene <- ClosestFeature(
  object = l1_teaseq_obj,
  regions    = dap_gr,
  annotation = Annotation(l1_teaseq_obj)
)
l1_teaseq_unimodal_atac_clust_dap <- l1_teaseq_unimodal_atac_clust_dap %>%
  mutate(
    closest_gene_symbol = dap_closest_gene$gene_name,
    closest_ensembl_id = dap_closest_gene$gene_id
  )
l1_teaseq_unimodal_atac_clust_dap_markers <- l1_teaseq_unimodal_atac_clust_dap %>%
  rename(
    l1_unimodal_atac_clust_numb = cluster,
    atac_peak = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_teaseq_unimodal_atac_clust_dap_markers <- l1_teaseq_unimodal_atac_clust_dap_markers %>% arrange(l1_unimodal_atac_clust_numb, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_teaseq_unimodal_atac_clust_dap_markers <- l1_teaseq_unimodal_atac_clust_dap_markers %>% relocate(c('l1_unimodal_atac_clust_numb','l1_unimodal_atac_clust_annot','atac_peak','closest_gene_symbol','closest_ensembl_id','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_teaseq_unimodal_atac_clust_dap_markers, '/filepath/fig1/fig1_supp/l1_unimodal_clust_annotation/l1_teaseq_unimodal_atac_clust_dap_markers.rds')

# TEAseq Unimodal ATAC (GeneActivity) L1 Cluster Markers
Idents(l1_teaseq_obj) <- 'atac.bulk.cluster.harmony'
DefaultAssay(l1_teaseq_obj) <- 'ACT'
l1_teaseq_unimodal_atac_clust_dag <- FindAllMarkers(l1_teaseq_obj, assay = 'ACT') # default parameters
l1_teaseq_unimodal_atac_clust_dag$l1_unimodal_atac_clust_annot <- l1_atac_clust_names[as.character(l1_teaseq_unimodal_atac_clust_dag$cluster)]
l1_teaseq_unimodal_atac_clust_dag_markers <- l1_teaseq_unimodal_atac_clust_dag %>%
  rename(
    l1_unimodal_atac_clust_numb = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_teaseq_unimodal_atac_clust_dag_markers <- l1_teaseq_unimodal_atac_clust_dag_markers %>% arrange(l1_unimodal_atac_clust_numb, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_teaseq_unimodal_atac_clust_dag_markers <- l1_teaseq_unimodal_atac_clust_dag_markers %>% relocate(c('l1_unimodal_atac_clust_numb','l1_unimodal_atac_clust_annot','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_teaseq_unimodal_atac_clust_dag_markers, '/filepath/fig1/fig1_supp/l1_unimodal_clust_annotation/l1_teaseq_unimodal_atac_clust_dag_markers.rds')

# TEAseq Unimodal RNA L1 Cluster Markers
DefaultAssay(l1_teaseq_obj) <- 'RNA'
l1_teaseq_unimodal_rna_clust_deg <- FindAllMarkers(l1_teaseq_obj, assay = 'RNA') # default parameters
l1_teaseq_unimodal_rna_clust_deg$l1_unimodal_rna_clust_annot <- l1_rna_clust_names[as.character(l1_teaseq_unimodal_rna_clust_deg$cluster)]
l1_teaseq_unimodal_rna_clust_deg_markers <- l1_teaseq_unimodal_rna_clust_deg %>%
  rename(
    l1_unimodal_rna_clust_numb = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_teaseq_unimodal_rna_clust_deg_markers <- l1_teaseq_unimodal_rna_clust_deg_markers %>% arrange(l1_unimodal_rna_clust_numb, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_teaseq_unimodal_rna_clust_deg_markers <- l1_teaseq_unimodal_rna_clust_deg_markers %>% relocate(c('l1_unimodal_rna_clust_numb','l1_unimodal_rna_clust_annot','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_teaseq_unimodal_rna_clust_deg_markers, '/filepath/fig1/fig1_supp/l1_unimodal_clust_annotation/l1_teaseq_unimodal_rna_clust_deg_markers.rds')

# TEAseq Unimodal ADT Cluster Markers - ADT (CLR-normalized)
DefaultAssay(l1_teaseq_obj) <- 'ADT'
Idents(l1_teaseq_obj) <- 'adt.bulk.cluster.harmony'
l1_teaseq_unimodal_adt_clust_dep <- FindAllMarkers(l1_teaseq_obj, assay = 'ADT', logfc.threshold = 0, min.pct = 0)
l1_teaseq_unimodal_adt_clust_dep$l1_unimodal_adt_clust_annot <- l1_adt_clust_names[as.character(l1_teaseq_unimodal_adt_clust_dep$cluster)]
l1_teaseq_unimodal_adt_clust_dep_markers <- l1_teaseq_unimodal_adt_clust_dep %>%
  rename(
    l1_unimodal_adt_clust_numb = cluster,
    adt_epitope = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_teaseq_unimodal_adt_clust_dep_markers <- l1_teaseq_unimodal_adt_clust_dep_markers %>% arrange(l1_unimodal_adt_clust_numb, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_teaseq_unimodal_adt_clust_dep_markers$adt_epitope <- gsub("^adt-", "", l1_teaseq_unimodal_adt_clust_dep_markers$adt_epitope)
l1_teaseq_unimodal_adt_clust_dep_markers <- l1_teaseq_unimodal_adt_clust_dep_markers %>% relocate(c('l1_unimodal_adt_clust_numb','l1_unimodal_adt_clust_annot','adt_epitope','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_teaseq_unimodal_adt_clust_dep_markers, '/filepath/fig1/fig1_supp/l1_unimodal_clust_annotation/l1_teaseq_unimodal_adt_clust_dep_markers.rds')

# Export all cluster marker files as sheets within one XLSX spreadsheet
clust_marker_df_dir   <- '/filepath/fig1/fig1_supp/l1_unimodal_clust_annotation/'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) |> str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/fig1/fig1_supp/l1_unimodal_clust_annotation/teaseq_l1_unimodal_cluster_markers_per_modality.xlsx")  

# UMAP embeddings for L1 unimodal ATAC vs RNA vs ADT clustering, colored by Tonsil vs PBMC sample origin ----

# L1 ATAC (peak-based) UMAP
l1_atac_umap_tissue <- DimPlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = FALSE, label.box = TRUE, label.size = 6, pt.size = 0.5, repel = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_atac_umap_tissue
ggsave("/filepath/fig1/fig1_supp/l1_umap_color_by_tissue/l1_atac_umap_tissue.pdf", plot = l1_atac_umap_tissue, width = 6, height = 6)

# L1 RNA UMAP, rotated as noted above
l1_rna_umap_tissue_rot <- DimPlot(l1_teaseq_obj, reduction = "umap.rna.harmony_rot", label = FALSE, label.box = TRUE, label.size = 6, pt.size = 0.5, repel = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_rna_umap_tissue_rot
ggsave("/filepath/fig1/fig1_supp/l1_umap_color_by_tissue/l1_rna_umap_tissue_rot.pdf", plot = l1_rna_umap_tissue_rot, width = 6, height = 6)

# L1 ADT UMAP
l1_adt_umap_tissue <- DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = FALSE, label.box = TRUE, label.size = 6, pt.size = 0.5, repel = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_adt_umap_tissue
ggsave("/filepath/fig1/fig1_supp/l1_umap_color_by_tissue/l1_adt_umap_tissue.pdf", plot = l1_adt_umap_tissue, width = 6, height = 6)