# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S4
# Fig. S4 - Gene expression modality weights vary between cells of each TEAseq 3WNN L1 cluster. 

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

# Making Ternary plots for WNN modality weights for all cells, split by cluster (related to Fig. 1E, which only visualized cluster centroids) ----

# Get metadata from object, including modality weights output by the Seurat FindMultimodalNeighbors function
l1_teaseq_obj.md <- l1_teaseq_obj@meta.data
l1_teaseq_obj.md <- l1_teaseq_obj.md %>%
  rename(RNA = RNA.weight, ATAC = ATAC.weight, ADT = ADT.weight)
levels(l1_teaseq_obj.md$l1_wnn_annot) <- c("CD4 Tn","CD4 Tn Ribo","CD4 Tm","Treg","Tfh-like","CD8 Tn","CD8 Tm/Inn","CD8 UTC","NK","Myl","Myl CD16","NBC","MBC", "GCB","ASC")

# Set colors for visualization
l1_clust_cols_ordered <- setNames(
  c('#8ab1cc','#e0c1b4','#94a890','#c9a3c1','#e49a78',
    '#c49980','#add5df','#d699ac','#728782','#b6b7ba',
    '#e1d4c7','#657ab0','#d4b5e4','#d66d97','#dfcc78'),
  c("CD4 Tn","CD4 Tn Ribo","CD4 Tm","Treg","Tfh-like",
    "CD8 Tn","CD8 Tm/Inn","CD8 UTC","NK","Myl",
    "Myl CD16","NBC","MBC","GCB","ASC")
)

# Create ternary plots
l1_clust_mod_weights_tern_all_cells <- ggtern(data = l1_teaseq_obj.md, aes(x = ATAC, y = ADT, z = RNA)) +
  scale_fill_manual(
    name   = "Cluster",
    values = l1_clust_cols_ordered) +
  geom_point(aes(fill = l1_wnn_annot), size = 2, alpha = 1, shape = 21, color = 'black', stroke = 0.1) +
  theme_hidegrid_minor() +
  theme_bw(base_family = 'sans') +
  theme_showarrows() +
  theme_arrowlarge() +
  geom_confidence_tern(color = 'black') +
  facet_wrap2(~ l1_wnn_annot, ncol = 5,
              strip = strip_themed(
                background_x = elem_list_rect(fill = l1_clust_cols_ordered,
                                              color = rep("black", length(l1_clust_cols_ordered)),
                                              size  = rep(1, length(l1_clust_cols_ordered))))) + 
  theme_hidetitles() +
  theme(
    text               = element_text(family = "sans"),
    axis.title         = element_text(family = "sans", color = 'black'),
    axis.text          = element_text(family = "sans", color = 'black'),
    legend.text        = element_text(family = "sans", color = 'black'),
    legend.title       = element_text(family = "sans", color = 'black'),
    strip.text         = element_text(family = "sans", size = 20, color = 'black'),
    tern.axis.text =  element_text(family = 'sans', size= 12),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "black"),
    panel.grid.minor = element_line(color = "black"),
    tern.axis.arrow = element_line(linewidth = 3, color = 'black'),
    tern.axis.arrow.sep  = 0.175,
    tern.axis.vshift     = 0.1,
    panel.spacing = unit(0.25, "cm"),
    panel.spacing.y = unit(0.5, "cm"),
    tern.axis.text = element_text(family = 'sans', size= 12)) + NoLegend()

# Save plots
pdf("/filepath/fig1/fig1_supp/l1_wnn_tern_plots/l1_clust_mod_weights_tern_all_cells.pdf", width = 12, height = 8)
l1_clust_mod_weights_tern_all_cells
dev.off()