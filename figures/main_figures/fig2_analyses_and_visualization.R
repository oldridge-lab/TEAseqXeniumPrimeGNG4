# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. 2
# Fig. 2 - GC versus nonGC-like Tfh states are skewed in helper-polarity phenotype.

# Set up R working environment ----

# Set working directory to Fig 2
setwd('/filepath/fig2/')

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
library(writexl) # 1.5.4
library(ggnewscale) # 0.5.0
library(readxl) # 1.4.3
library(UCell) # 2.8.0
library(ggrepel) # 0.9.5

# Import Seurat object from Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors)
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Set Tfh subcluster annotations ----

# Mapping of 3WNN cluster numbers to cluster annotations
tfh_subclust_names <- c(
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

# Specify Seurat object cell idents
l3_teaseq_tfh_obj$tfh_wnn_annot <- l3_teaseq_tfh_obj$wnn.tfh.subcluster.harmony
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
l3_teaseq_tfh_obj <- RenameIdents(l3_teaseq_tfh_obj, tfh_subclust_names)
l3_teaseq_tfh_obj$tfh_wnn_annot <- Idents(l3_teaseq_tfh_obj)

# Set ident order for visualization
wnn_tfh_ident_order <- c('Tcm','Tfh-AP1','Tfh-Circ','Tfh-IL10','Tfh-Resting','Tfh-Int','Tfh-NFATC1','Tfh-CXCL13','Tfh-BOB1')
Idents(l3_teaseq_tfh_obj) <- factor(
  Idents(l3_teaseq_tfh_obj),
  levels = wnn_tfh_ident_order
)

# Annotating GC vs nonGC-like states
gc_vs_nongc_like <- c(
  "Tfh-Circ" = 'nonGC',
  "Tcm" = "nonGC",
  "Tfh-Int" = "GC",
  "Tfh-NFATC1" = "GC",
  "Tfh-IL10" = "GC",
  "Tfh-Resting" = "nonGC",
  "Tfh-CXCL13" = "GC",
  "Tfh-BOB1" = "GC",
  "Tfh-AP1" = "nonGC"
)
l3_teaseq_tfh_obj$gc_vs_nongc_like <- l3_teaseq_tfh_obj$tfh_wnn_annot
Idents(l3_teaseq_tfh_obj) <- 'gc_vs_nongc_like'
l3_teaseq_tfh_obj <- RenameIdents(l3_teaseq_tfh_obj, gc_vs_nongc_like)
l3_teaseq_tfh_obj$gc_vs_nongc_like <- Idents(l3_teaseq_tfh_obj)
table(l3_teaseq_tfh_obj$gc_vs_nongc_like)

# Specify colors for L3 Tfh-like clusters
tfh_clust_cols <- c(
  "Tfh-Circ" = '#657ab0',
  "Tcm" = "#728782",
  "Tfh-Int" = "#c49980",
  "Tfh-NFATC1" = "#c9a3c1",
  "Tfh-IL10" = "#dfcc78",
  "Tfh-Resting" = "#b6b7ba",
  "Tfh-CXCL13" = "#bad4f4",
  "Tfh-BOB1" = "#c9744d",
  "Tfh-AP1" = "#e1d4c7"
)

# Specify colors for each tissue
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# Fig 2A - Experimental schematic for parallel TEAseq and spectral flow cytometry analysis of Tfh cell states ----

# Created in BioRender and exported as PDF, edited it Illustrator

# Fig 2B - 3WNN Clustering and UMAP of L3 Tfh-like cells ----

tfh_wnn_umap_clust_dimplot <- DimPlot(l3_teaseq_tfh_obj, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 5.25, reduction = "umap.wnn.harmony", pt.size = 2, group.by = 'tfh_wnn_annot', cols = tfh_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig2/fig2b/tfh_wnn_umap_clust_dimplot.pdf', width = 5.5, height = 4)
tfh_wnn_umap_clust_dimplot
dev.off()

# Fig 2C - Barplot of L3 cluster proportion for each tissue type ----
tfh.md <- l3_teaseq_tfh_obj@meta.data
tissue_cell_tot <- tfh.md %>%
  group_by(hto.tissue) %>%
  summarise(total = n(), .groups = "drop")
tissue_cell_prop_df <- tfh.md %>%
  group_by(tfh_wnn_annot, hto.tissue) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(tissue_cell_tot, by = "hto.tissue") %>%
  mutate(prop = count / total)
tfh_clust_prop_per_tissue_barplot <- ggplot(tissue_cell_prop_df, aes(x = hto.tissue, y = prop, fill = tfh_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = tfh_clust_cols) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  labs(y = "Cluster Proportion") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 30, 'sans', margin = margin(t = 10)),
    axis.text.y = element_text(size = 20, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  ) +
  NoLegend()
tfh_clust_prop_per_tissue_barplot
ggsave("/filepath/fig2/fig2c/tfh_clust_prop_per_tissue_barplot.pdf", plot = tfh_clust_prop_per_tissue_barplot, width = 4.5, height = 6)

# Fig 2D - 3WNN L3 UMAP split by tissue ----

tfh_wnn_umap_tissue_dimplot <- DimPlot(l3_teaseq_tfh_obj, label = FALSE, reduction = "umap.wnn.harmony", pt.size = 2, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig2/fig2d/tfh_wnn_umap_tissue_dimplot.pdf', width = 5.5, height = 4)
tfh_wnn_umap_tissue_dimplot
dev.off()

# Fig 2E - UMAP visualization of differentially expressed features across L3 Tfh-like cell subclusters from ADT, RNA, and SCENIC (RNA-based) assays ----

# ADT Feature Plots

DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
pdf('/filepath/fig2/fig2e/adt/tfh_cd127_adt_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD127',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
pdf('/filepath/fig2/fig2e/adt/tfh_cd25_adt_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD25',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
pdf('/filepath/fig2/fig2e/adt/tfh_cd69_adt_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD69',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
pdf('/filepath/fig2/fig2e/adt/tfh_pd1_adt_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD279',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
pdf('/filepath/fig2/fig2e/adt/tfh_tigit_adt_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'adt-TIGIT',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# RNA Feature Plots

DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig2/fig2e/rna/tfh_tox2_rna_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'TOX2',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig2/fig2e/rna/tfh_il10_rna_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'IL10',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig2/fig2e/rna/tfh_il21_rna_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'IL21',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig2/fig2e/rna/tfh_il4_rna_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'IL4',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig2/fig2e/rna/tfh_cxcl13_rna_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'CXCL13',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# SCENIC Feature Plots

DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
pdf('/filepath/fig2/fig2e/scenic/tfh_klf2_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'KLF2-extended (61g)',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
pdf('/filepath/fig2/fig2e/scenic/tfh_prdm1_scenic_umap.pdf', width = 6, height = 4) # BLIMP1
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'PRDM1-extended (32g)',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
pdf('/filepath/fig2/fig2e/scenic/tfh_ikaros_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'IKZF1-extended (317g)',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
pdf('/filepath/fig2/fig2e/scenic/tfh_nfatc1_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'NFATC1-extended (184g)',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
pdf('/filepath/fig2/fig2e/scenic/tfh_bcl6_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'BCL6-extended (105g)',
               use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, 
               legend.position = 'right', legend.length = 11, border.size = 1, pt.size = 2) +
  theme(legend.position = c(1, 0.5),
        legend.direction = 'vertical',
        plot.margin = margin(r = 10)) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig 2F - Dot plot visualizing T helper features across Tfh-like cells only ----

# Visualizing Tfh clusters only - excluding Tcm cluster by filtering L3 object and creating new L4 Tfh-only object
l4_teaseq_tfh_only_obj <- subset(l3_teaseq_tfh_obj, tfh_wnn_annot != 'Tcm')

# Calculate Th17 RNA gene signature enrichment using UCell
th17_gene_sig <- c('RORA','RORC','CCR6','KLRB1','IL23R','IL17A','IL17F','IL22','IL26','CCL20','RUNX1','PALLD','TNFRSF8') # gene signature curated from human Th17 literature review, refer to Main Text and Supplementary Materials for references
th17_gene_sig_list <- list(
  Th17 = th17_gene_sig
)
DefaultAssay(l4_teaseq_tfh_only_obj) <- 'RNA'
l4_teaseq_tfh_only_obj <- AddModuleScore_UCell(l4_teaseq_tfh_only_obj, features= th17_gene_sig_list, name= NULL)

# Create dot plot of polarity ADT markers and Th17 RNA gene signature
DefaultAssay(l4_teaseq_tfh_only_obj) <- 'ADT'
Idents(l4_teaseq_tfh_only_obj) <- 'tfh_wnn_annot'
wnn_tfh_ident_order <- c('Tfh-BOB1','Tfh-CXCL13','Tfh-NFATC1','Tfh-Int','Tfh-IL10','Tfh-AP1','Tfh-Resting','Tfh-Circ')
Idents(l4_teaseq_tfh_only_obj) <- factor(
  Idents(l4_teaseq_tfh_only_obj),
  levels = wnn_tfh_ident_order
)
tfh_polarity_dotplot_feats <- c('adt-CD183','adt-CD194','adt-CD196','adt-CD161','Th17')
tfh_polarity_dotplot <- DotPlot(l4_teaseq_tfh_only_obj, features = tfh_polarity_dotplot_feats,
                                              cols = 'PuOr', scale.by = 'size', dot.scale = 10) + 
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10, 'sans'),
    axis.text.y      = element_text(size = 10, 'sans'),
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    legend.position  = 'bottom',
    legend.justification = c(1,0),
    legend.box       = 'horizontal',
    legend.key.size     = unit(0.1, "cm"),
    legend.text         = element_text(size = 10),
    legend.title        = element_text(size = 7),
    legend.spacing.x    = unit(0.1, "cm")
  ) +
  geom_vline(xintercept = seq(1.5, length(tfh_polarity_dotplot_feats) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(l4_teaseq_tfh_only_obj$tfh_wnn_annot))    - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  guides(
    size = guide_legend(
      title = "", # labels added manually in Illustrator
      order = 1,
      barwidth = 0.1
    ),
    color = guide_colorbar(
      title = "", # labels added manually in Illustrator
      order = 2,
      barwidth = 5,
      barheight = 1
    )
  )  + 
  scale_x_discrete(labels = c('CXCR3 [ADT]','CCR4 [ADT]','CCR6 [ADT]','KLRB1 [ADT]','Th17 Sig [RNA]'))
pdf('/filepath/fig2/fig2f/tfh_polarity_dotplot.pdf', height = 4.5, width = 3.75)
tfh_polarity_dotplot
dev.off()

# Fig 2G - Contour plots of spectral flow cytometry gating scheme for Tfh subsets in Tonsil versus PBMC ----

# Created in FlowJo, annotated in Illustrator

# Fig 2H - Barplots of TF GMFI in each gated CD4 T cell subset ----

# Fig 2H - PBMC T-bet GMFI Barplot
tbet_gmfi_tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "Tbet")
colnames(tbet_gmfi_tfh_flow_df) <- c('Tissue','Donor','Naïve','Th1','cTfh0','cTfh1','cTfh2','cTfh17')
tfh_barplot_cols <- c(
  "Naïve" = 'grey',
  "cTfh0" = '#7570b3',
  "cTfh1" = "#ff7e79",
  "cTfh2" = "#009091",
  "cTfh17" = "#fbb041",
  "Th1" = "coral2",
  "Th2" = "#639091",
  "Th17" = "orange3"
)
tbet_gmfi_tfh_flow_df <- tbet_gmfi_tfh_flow_df %>% subset(Tissue == 'PBMC')
tbet_gmfi_tfh_flow_df_long <- tbet_gmfi_tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:cTfh17,               
    names_to  = "Subset",                  
    values_to = "gMFI"                     
  )
tbet_gmfi_tfh_flow_df_long <- tbet_gmfi_tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c('Naïve','Th1','cTfh0','cTfh1','cTfh2','cTfh17')
    )
  )
tbet_gmfi_tfh_barplot <- ggplot(tbet_gmfi_tfh_flow_df_long, aes(x = Subset, y = gMFI, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal, 
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = tbet_gmfi_tfh_flow_df_long,
    shape     = 21,
    color   = "black",     # outline color
    stroke  = 0.5,         # thin border
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tfh_barplot_cols) +
  theme(axis.title = element_blank())
tbet_gmfi_tfh_barplot
pdf('/filepath/fig2/fig2h/pbmc_tbet_gmfi_tfh_barplot.pdf', width = 5, height = 5)
tbet_gmfi_tfh_barplot
dev.off()

# Fig 2H - Tonsil T-bet GMFI Barplot
tbet_gmfi_tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "Tbet")
tfh_barplot_cols <- c(
  "Naïve" = 'grey',
  "Tfh0" = '#7570b3',
  "Tfh1" = "#ff7e79",
  "Tfh2" = "#009091",
  "Tfh17" = "#fbb041",
  "Th1" = "coral2",
  "Th2" = "#639091",
  "Th17" = "orange3"
)
tbet_gmfi_tfh_flow_df <- tbet_gmfi_tfh_flow_df %>% subset(Tissue == 'Tonsil')
tbet_gmfi_tfh_flow_df_long <- tbet_gmfi_tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:Tfh17,               
    names_to  = "Subset",                  
    values_to = "gMFI"                     
  )
tbet_gmfi_tfh_flow_df_long <- tbet_gmfi_tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c("Naïve","Th1","Tfh0", "Tfh1", "Tfh2", "Tfh17")
    )
  )
tbet_gmfi_tfh_barplot <- ggplot(tbet_gmfi_tfh_flow_df_long, aes(x = Subset, y = gMFI, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal, 
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = tbet_gmfi_tfh_flow_df_long,
    shape     = 21,
    color   = "black",     # outline color
    stroke  = 0.5,         # thin border
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tfh_barplot_cols) +
  theme(axis.title = element_blank())
tbet_gmfi_tfh_barplot
pdf('/filepath/fig2/fig2h/tonsil_tbet_gmfi_tfh_barplot.pdf', width = 5, height = 5)
tbet_gmfi_tfh_barplot
dev.off()

# Fig 2H - PBMC GATA3 GMFI Barplot
gata3_gmfi_tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "GATA3")
colnames(gata3_gmfi_tfh_flow_df) <- c('Tissue','Donor','Naïve','Th2','cTfh0','cTfh1','cTfh2','cTfh17')
tfh_barplot_cols <- c(
  "Naïve" = 'grey',
  "cTfh0" = '#7570b3',
  "cTfh1" = "#ff7e79",
  "cTfh2" = "#009091",
  "cTfh17" = "#fbb041",
  "Th1" = "coral2",
  "Th2" = "#639091",
  "Th17" = "orange3"
)
gata3_gmfi_tfh_flow_df <- gata3_gmfi_tfh_flow_df %>% subset(Tissue == 'PBMC')
gata3_gmfi_tfh_flow_df_long <- gata3_gmfi_tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:cTfh17,               
    names_to  = "Subset",                  
    values_to = "gMFI"                     
  )
gata3_gmfi_tfh_flow_df_long <- gata3_gmfi_tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c('Naïve','Th2','cTfh0','cTfh1','cTfh2','cTfh17')
    )
  )
gata3_gmfi_tfh_barplot <- ggplot(gata3_gmfi_tfh_flow_df_long, aes(x = Subset, y = gMFI, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal, 
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = gata3_gmfi_tfh_flow_df_long,
    shape     = 21,
    color   = "black",     # outline color
    stroke  = 0.5,         # thin border
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tfh_barplot_cols) +
  theme(axis.title = element_blank())
gata3_gmfi_tfh_barplot
pdf('/filepath/fig2/fig2h/pbmc_gata3_gmfi_tfh_barplot.pdf', width = 5, height = 5)
gata3_gmfi_tfh_barplot
dev.off()

# Fig 2H - Tonsil GATA3 GMFI Barplot
gata3_gmfi_tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "GATA3")
tfh_barplot_cols <- c(
  "Naïve" = 'grey',
  "Tfh0" = '#7570b3',
  "Tfh1" = "#ff7e79",
  "Tfh2" = "#009091",
  "Tfh17" = "#fbb041",
  "Th1" = "coral2",
  "Th2" = "#639091",
  "Th17" = "orange3"
)
gata3_gmfi_tfh_flow_df <- gata3_gmfi_tfh_flow_df %>% subset(Tissue == 'Tonsil')
gata3_gmfi_tfh_flow_df_long <- gata3_gmfi_tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:Tfh17,               
    names_to  = "Subset",                  
    values_to = "gMFI"                     
  )
gata3_gmfi_tfh_flow_df_long <- gata3_gmfi_tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c("Naïve","Th2","Tfh0", "Tfh1", "Tfh2", "Tfh17")
    )
  )
gata3_gmfi_tfh_barplot <- ggplot(gata3_gmfi_tfh_flow_df_long, aes(x = Subset, y = gMFI, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal, 
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = gata3_gmfi_tfh_flow_df_long,
    shape     = 21,
    color   = "black",     # outline color
    stroke  = 0.5,         # thin border
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tfh_barplot_cols) +
  theme(axis.title = element_blank())
gata3_gmfi_tfh_barplot
pdf('/filepath/fig2/fig2h/tonsil_gata3_gmfi_tfh_barplot.pdf', width = 5, height = 5)
gata3_gmfi_tfh_barplot
dev.off()

# Fig 2H - PBMC RORgt GMFI Barplot
rorgt_gmfi_tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "RORgt")
colnames(rorgt_gmfi_tfh_flow_df) <- c('Tissue','Donor','Naïve','Th17','cTfh0','cTfh1','cTfh2','cTfh17')
tfh_barplot_cols <- c(
  "Naïve" = 'grey',
  "cTfh0" = '#7570b3',
  "cTfh1" = "#ff7e79",
  "cTfh2" = "#009091",
  "cTfh17" = "#fbb041",
  "Th1" = "coral2",
  "Th2" = "#639091",
  "Th17" = "orange3"
)
rorgt_gmfi_tfh_flow_df <- rorgt_gmfi_tfh_flow_df %>% subset(Tissue == 'PBMC')
rorgt_gmfi_tfh_flow_df_long <- rorgt_gmfi_tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:cTfh17,               
    names_to  = "Subset",                  
    values_to = "gMFI"                     
  )
rorgt_gmfi_tfh_flow_df_long <- rorgt_gmfi_tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c('Naïve','Th17','cTfh0','cTfh1','cTfh2','cTfh17')
    )
  )
rorgt_gmfi_tfh_barplot <- ggplot(rorgt_gmfi_tfh_flow_df_long, aes(x = Subset, y = gMFI, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal, 
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = rorgt_gmfi_tfh_flow_df_long,
    shape     = 21,
    color   = "black",     # outline color
    stroke  = 0.5,         # thin border
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tfh_barplot_cols) +
  theme(axis.title = element_blank())
rorgt_gmfi_tfh_barplot
pdf('/filepath/fig2/fig2h/pbmc_rorgt_gmfi_tfh_barplot.pdf', width = 5, height = 5)
rorgt_gmfi_tfh_barplot
dev.off()

# Fig 2H - Tonsil RORgt GMFI Barplot
rorgt_gmfi_tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "RORgt")
tfh_barplot_cols <- c(
  "Naïve" = 'grey',
  "Tfh0" = '#7570b3',
  "Tfh1" = "#ff7e79",
  "Tfh2" = "#009091",
  "Tfh17" = "#fbb041",
  "Th1" = "coral2",
  "Th2" = "#639091",
  "Th17" = "orange3"
)
rorgt_gmfi_tfh_flow_df <- rorgt_gmfi_tfh_flow_df %>% subset(Tissue == 'Tonsil')
rorgt_gmfi_tfh_flow_df_long <- rorgt_gmfi_tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:Tfh17,               
    names_to  = "Subset",                  
    values_to = "gMFI"                     
  )
rorgt_gmfi_tfh_flow_df_long <- rorgt_gmfi_tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c("Naïve","Th17","Tfh0", "Tfh1", "Tfh2", "Tfh17")
    )
  )
rorgt_gmfi_tfh_barplot <- ggplot(rorgt_gmfi_tfh_flow_df_long, aes(x = Subset, y = gMFI, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal, 
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = rorgt_gmfi_tfh_flow_df_long,
    shape     = 21,
    color   = "black",     # outline color
    stroke  = 0.5,         # thin border
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tfh_barplot_cols) +
  theme(axis.title = element_blank())
rorgt_gmfi_tfh_barplot
pdf('/filepath/fig2/fig2h/tonsil_rorgt_gmfi_tfh_barplot.pdf', width = 5, height = 5)
rorgt_gmfi_tfh_barplot
dev.off()

# Fig 2I - Contour plots for spectral flow cytometry gating of TIGIT vs CD127 in PBMC cTfh helper-polarized subsets ----

# PDF exported from FlowJo, further annotated in Illustrator

# Fig 2J - Contour plots for spectral flow cytometry gating of TIGIT vs CD127 in Tonsil Tfh helper-polarized subsets ----

# PDF exported from FlowJo, further annotated in Illustrator

# Fig 2K - Stacked barplots of TIGIT vs CD127 quadrant subset frequency in helper-polarized Tfh and cTfh subsets ----

# Prepare dataframe for plotting
tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "TIGIT_IL7R_Pct")
tfh_flow_df_long <- tfh_flow_df %>%
  pivot_longer(
    cols       = -Donor,                       
    names_to   = c("Subset", "Quartile"),      
    names_sep  = "_",
    values_to  = "Frequency"
  ) %>%
  mutate(
    Subset   = factor(Subset,   levels = c("Tfh0", "Tfh1", "Tfh2", "Tfh17","cTfh0", "cTfh1", "cTfh2", "cTfh17")),
    Quartile = factor(Quartile, levels = c("Q1", "Q2", "Q3", "Q4"))
  )
tfh_flow_df_long_summary <- tfh_flow_df_long %>%
  group_by(Subset, Quartile) %>%
  summarise(
    mean_freq = mean(Frequency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup()

# Specify quadrant colors
quadrant_pct_colors <- c(
  "Q1" = "#f09292",   # TIGIT+
  "Q2" = "#c9a3c1",   # DP
  "Q3" = "#bad4f4",   # IL7R+
  "Q4" = "#b6b7ba"    # DN
)

# PBMC cTfh subsets
pbmc_ctfh_cd127_tigit_summary_df <- tfh_flow_df_long_summary %>% subset(Subset %in% c('cTfh0','cTfh1','cTfh2','cTfh17'))
tfh_quadrant_freq_barplot <- ggplot(pbmc_ctfh_cd127_tigit_summary_df, aes(x = Subset, y = mean_freq, fill = Quartile)) +
  geom_bar(
    stat  = "identity",
    color = "black",    
    width = 0.8
  ) +
  scale_fill_manual(values = quadrant_pct_colors, name = "Quartile") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    axis.text.x     = element_blank(), # element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 18),
    axis.title.y     = element_blank(),  # element_text(size = 20),
    axis.title.x      = element_blank(),
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 11)
  ) +
  labs(
    y     = "Quadrant Frequency"
  ) + NoLegend()
tfh_quadrant_freq_barplot
pdf('/filepath/fig2/fig2k/pbmc_ctfh_cd127_tigit_quadrant_freq_barplot.pdf', width = 4, height = 4)
tfh_quadrant_freq_barplot
dev.off()

# Tonsil Tfh subsets
tonsil_tfh_cd127_tigit_summary_df <- tfh_flow_df_long_summary %>% subset(Subset %in% c('Tfh0','Tfh1','Tfh2','Tfh17'))
tfh_quadrant_freq_barplot <- ggplot(tonsil_tfh_cd127_tigit_summary_df, aes(x = Subset, y = mean_freq, fill = Quartile)) +
  geom_bar(
    stat  = "identity",
    color = "black",    
    width = 0.8
  ) +
  scale_fill_manual(values = quadrant_pct_colors, name = "Quartile") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    axis.text.x     = element_blank(), # element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 18),
    axis.title.y     = element_blank(),  # element_text(size = 20),
    axis.title.x      = element_blank(),
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 11)
  ) +
  labs(
    #  x     = "T Cell Subset",
    y     = "Quadrant Frequency",
  ) + NoLegend()
tfh_quadrant_freq_barplot
pdf('/filepath/fig2/fig2k/tonsil_tfh_cd127_tigit_quadrant_freq_barplot.pdf', width = 4, height = 4)
tfh_quadrant_freq_barplot
dev.off()

# Fig 2L - Contour plots of chemokine receptor versus BCL6 expression in polarized Tfh ----

# PDF exported from FlowJo, further annotated in Illustrator

# Fig 2M -  Barplot of percent Bcl6 positive across CD4 T cell subsets in tonsil ----

tfh_flow_df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "Bcl6_Pct")
tfh_barplot_cols <- c(
  "Naïve"  = "grey",
  "nonTfh" = "#c49980",
  "Tfh0"    = "#7570b3",
  "Tfh1"    = "#ff7e79",
  "Tfh2"    = "#009091",
  "All"   = "#fbb040",
  'IL7R+ TIGIT-' = '#f2d48d',
  'IL7R- TIGIT+' = 'darkorange3'
)
tfh_flow_df_long <- tfh_flow_df %>%
  pivot_longer(
    cols      = Naïve:`IL7R- TIGIT+`,               
    names_to  = "Subset",                  
    values_to = "Freq"                     
  )
tfh_flow_df_long <- tfh_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c("Naïve", "nonTfh", "Tfh0", "Tfh1", "Tfh2", "All","IL7R+ TIGIT-", 'IL7R- TIGIT+')
    )
  )
tfh_bcl6_pct_barplot <- ggplot(tfh_flow_df_long, aes(x = Subset, y = Freq, fill = Subset)) +
  stat_summary(
    fun    = mean,
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) + 
  geom_jitter(
    data      = tfh_flow_df_long,
    shape     = 21,
    color   = "black",    
    stroke  = 0.5,       
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.8,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  labs(
    y     = "% BCL6+ of Subset") +
  scale_fill_manual(values = tfh_barplot_cols) +
  scale_x_discrete(labels = c("Naïve", "nonTfh", "Tfh0", "Tfh1", "Tfh2", "","",""))
tfh_bcl6_pct_barplot
pdf('/filepath/fig2/fig2l/tfh_bcl6_pct_barplot.pdf', width = 6.5, height = 5)
tfh_bcl6_pct_barplot
dev.off()

# Data file S10 - 95% confidence interval calculations for TF GMFI and BCL6 expression percentage analysis of spectral flow cytometry data ----

# Select sheets from raw flow data spreadsheet
xlsx_path <- "/filepath/fig2/fig2_flow_data_spreadsheet.xlsx"
sheets <- c("Tbet", "GATA3", "RORgt", "Bcl6_Pct")

# Function to retrieve tissue and donor groups from raw data, then compute relevant 95CI metrics
process_sheet <- function(sheet_name) {
  
  df <- read_excel(xlsx_path, sheet = sheet_name)
  subset_cols <- df %>%
    select(-any_of(c("Tissue","Donor"))) %>%
    select(where(is.numeric)) %>%
    names()
  
  df %>%
    pivot_longer(all_of(subset_cols), names_to = "Subset", values_to = "Value") %>%
    group_by(Tissue, Subset) %>%
    summarise(
      n = n(),
      mean = mean(Value),
      sd = sd(Value),
      se = sd / sqrt(n),
      ci = qt(0.975, df = n - 1) * se,
      ymin = mean - ci,
      ymax = mean + ci,
      .groups = "drop"
    ) %>%
    mutate(sheet = sheet_name, .before = 1)
}

# Execute function
all_stats <- map_dfr(sheets, process_sheet)

# Save 95CI metrics
writexl::write_xlsx(
  list(Stats_by_Group = all_stats %>%
         select(sheet, Tissue, Subset, n, mean, sd, se, ci, ymin, ymax)),
  path = "/filepath/fig2/fig2_flow_data_95ci_stats.xlsx"
)

# Data file S10 - 95% confidence interval calculations for TIGIT vs CD127 quadrant percentage analysis of spectral flow cytometry data ----

# Prepare dataframe
df <- read_excel("/filepath/fig2/fig2_flow_data_spreadsheet.xlsx", sheet = "TIGIT_IL7R_Pct")
df_long <- df %>%
  pivot_longer(-Donor, names_to = "Var", values_to = "Percent") %>%
  separate(Var, into = c("Subset","Quadrant"), sep = "_Q") %>%
  mutate(Quadrant = paste0("Q", Quadrant))

# Compute mean ± 95% CI across donors
tigit_cd127_quadrants_stats <- df_long %>%
  group_by(Subset, Quadrant) %>%
  summarise(
    n = n(),
    mean = mean(Percent),
    sd = sd(Percent),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    ymin = mean - ci,
    ymax = mean + ci,
    .groups = "drop"
  )

# Save spreadsheet
writexl::write_xlsx(
  tigit_cd127_quadrants_stats,
  path = "/filepath/fig2/fig2_flow_data_tigit_cd127_quadrants.xlsx"
)