# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S7
# Fig. S7 - Tfh adopt distinct states in tonsil and peripheral blood unified by a set of multimodal gene expression features.

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

# Import Seurat object from TEAseq Data Preprocessing Step 14 (trimodal dimensionality reduction and L2 3WNN subclustering of T cells from L1 object, including Harmony integration across donors)
l2_teaseq_tcell_obj <- readRDS('/filepath/step14_tcell_subcluster/l2_teaseq_tcell_obj.rds')

# L2 T cell TEAseq object cluster names
l2_tcell_clust_names <- c(
  '0' = 'CD4 Tn 1',
  '1' = 'CD4 Tn 2',
  '2' = 'CD4 Tn 3',
  '3' = 'CD4 Tn 4',
  '4' = 'Tfh GC',
  '5' = 'CD4 Tcm/fh',
  '6' = 'CD8 Tn',
  '7' = 'CD4 Tem',
  '8' = 'CD8 CTL',
  '9' = 'CD4 Tn 5',
  '10' = 'cTreg',
  '11' = 'CD4 Tn 6',
  '12' = 'PLZF Inn',
  '13' = 'gdT',
  '14' = 'Th17',
  '15' = 'CD8 Tm',
  '16' = 'Treg RORgt',
  '17' = 'Tfh IL10',
  '18' = 'CD4 Tn 7',
  '19' = 'CD8 UTC',
  '20' = 'CD4 Inn'
)

# Label L2 3WNN clusters with assigned names
l2_teaseq_tcell_obj$l2_tcell_wnn_annot <- l2_teaseq_tcell_obj$wnn.tcell.subcluster.harmony
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'
l2_teaseq_tcell_obj <- RenameIdents(l2_teaseq_tcell_obj, l2_tcell_clust_names)
l2_teaseq_tcell_obj$l2_tcell_wnn_annot <- Idents(l2_teaseq_tcell_obj)
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'

# Specify L2 TEAseq object colors (T cell subclustering)
tcell_clust_cols <- c(
  "CD4 Tn 1"               = "#8ab1cc",
  "CD4 Tn 2"               = "#add5df",
  "CD4 Tn 3"               = "#94a890",
  "CD4 Tn 4"               = "#e0c1b4",
  "Tfh GC"                 = "#e49a78",
  "CD4 Tcm/fh"             = "#657ab0",
  "CD8 Tn"                 = "#e17794",
  "CD4 Tem"                = "#72d1b4",
  "CD8 CTL"                = "#c9744d",
  "CD4 Tn 5"               = "#d66d97",
  "cTreg"                  = "#b37cbf",
  "CD4 Tn 6"               = "#d4b5e4",
  "PLZF Inn"               = "#e1d4c7",
  "gdT"                    = "#c9a3c1",
  "Th17"                   = "#f4a6c7",
  "CD8 Tm"                 = "#bad4f4",
  "Treg RORgt"             = "#dfcc78",
  "Tfh IL10"               = "#728762",
  "CD4 Tn 7"               = "#c49980",
  "CD8 UTC"                = "#b6b7ba",
  "CD4 Inn"                = "#d699ac"
)

# A) L2 3WNN UMAP with T cell cluster labels ----
tcell_wnn_umap_coords <- FetchData(l2_teaseq_tcell_obj, vars = c("harmonywnnUMAP_1", "harmonywnnUMAP_2")) %>%
  mutate(l2_tcell_wnn_annot = Idents(l2_teaseq_tcell_obj),
         fill_col = tcell_clust_cols[as.character(l2_tcell_wnn_annot)])
l1_wnn_umap_lab_df <- tcell_wnn_umap_coords %>% 
  group_by(l2_tcell_wnn_annot, fill_col) %>% 
  summarize(UMAP_1 = median(harmonywnnUMAP_1), UMAP_2 = median(harmonywnnUMAP_2), .groups = 'drop')

Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'

tcell_wnn_umap_clust_lab <- DimPlot(l2_teaseq_tcell_obj, reduction = "umap.wnn.harmony", pt.size = 1.5, 
                                    group.by = 'l2_tcell_wnn_annot', cols = tcell_clust_cols) + 
  new_scale_color() + 
  geom_label_repel(
    data = l1_wnn_umap_lab_df,
    aes(x = UMAP_1, y = UMAP_2, label = l2_tcell_wnn_annot, fill = fill_col, color = 'black'),
    size = 6,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines"),
    force = 1,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  theme(plot.title = element_blank()) +
  scale_fill_identity() +
  scale_color_identity() +
  NoLegend() +
  NoAxes() +
  coord_fixed()
pdf('/filepath/fig1/fig1_supp/l2_tcell_wnn_umap_clust_lab.pdf', width = 8, height = 8)
tcell_wnn_umap_clust_lab
dev.off()

# B) L2 3WNN T Cell UMAP split by tissue ----
tcell_wnn_umap_tissue <- DimPlot(l2_teaseq_tcell_obj, reduction = "umap.wnn.harmony", label = FALSE, label.box = TRUE, label.size = 4, pt.size = 1.5, repel = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1/fig1_supp/tcell_wnn_umap_tissue.pdf', width = 8, height = 8)
tcell_wnn_umap_tissue
dev.off()

# C) Barplot of donor proportion per L2 3WNN cluster normalized to donor cell total ----

# Determine normalized and rescaled proportions per cluster only using unenriched samples (i.e. not CD4-enriched samples), as in Fig 1C
l2_teaseq_tcell_obj$cd4_vs_bulk <- ifelse(grepl("-CD4$", l2_teaseq_tcell_obj$hto.sort), "CD4-Enriched",
                                          ifelse(grepl("-Bulk$", l2_teaseq_tcell_obj$hto.sort), "All Cells", NA))
Idents(l2_teaseq_tcell_obj) <- 'cd4_vs_bulk'
tcell_unenr_samp_md <- subset(l2_teaseq_tcell_obj, idents = c('All Cells'))@meta.data # only retrieve samples containing all mononuclear cell types - not the '-CD4' suffix samples enriched for CD4+
donor_totals <- tcell_unenr_samp_md %>%
  group_by(hto.donor) %>%
  summarise(total = n(), .groups = "drop") # number of cells from each donor in each cluster
norm_prop_df <- tcell_unenr_samp_md %>%
  group_by(l2_tcell_wnn_annot, hto.donor) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(donor_totals, by = "hto.donor") %>%
  mutate(raw_fraction = count / total) # within donor proportion of cells corresponding to given cluster
norm_prop_df <- norm_prop_df %>%
  group_by(l2_tcell_wnn_annot) %>%
  mutate(sum_raw = sum(raw_fraction),
         normalized_fraction = raw_fraction / sum_raw) %>% # rescale proportions to sum to 1
  ungroup()

# Find clusters enriched in tonsil samples, and sort clusters from low to high for barplot visualization
t_proportions <- norm_prop_df %>%
  filter(hto.donor %in% paste0("T", 1:4)) %>%
  group_by(l2_tcell_wnn_annot) %>%
  summarise(t_sum = sum(normalized_fraction)) %>%
  complete(l2_tcell_wnn_annot = levels(norm_prop_df$l2_tcell_wnn_annot), fill = list(t_sum = 0)) %>%
  arrange(t_sum)
barplot_cluster_order <- t_proportions$l2_tcell_wnn_annot
norm_prop_df <- norm_prop_df %>% mutate(l2_tcell_wnn_annot = factor(l2_tcell_wnn_annot, levels = barplot_cluster_order))

# Specify colors for each tissue and donor
pbmc_colors <- colorRampPalette(c("steelblue4", "steelblue1"))(4)
tonsil_colors <- colorRampPalette(c("darkorange3", "orange"))(4)
names(pbmc_colors) <- paste0("P", 1:4)
names(tonsil_colors) <- paste0("T", 1:4)
donor_colors <- c(pbmc_colors, tonsil_colors)

# Assemble barplot
tcell_clust_prop_norm_barplot <- ggplot(norm_prop_df, aes(x = factor(l2_tcell_wnn_annot), y = normalized_fraction, fill = hto.donor)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = donor_colors) +
  labs(y = "Normalized Donor Proportion",
       fill = "Donor") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10), 'sans'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, angle = 45, 'sans', hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
tcell_clust_prop_norm_barplot
ggsave("/filepath/fig1/fig1_supp/tcell_clust_prop_norm_barplot.pdf", plot = tcell_clust_prop_norm_barplot, width = 12, height = 8)

# D) L2 3WNN UMAP T cell feature plots ----

DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cxcr5_prom_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'chr11-118883255-118884245', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_icos_enh_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'chr2-204022568-204023641', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cxcr5_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'CXCR5', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_tox2_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'TOX2', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_foxp3_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'FOXP3', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_klf2_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'KLF2', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cd69_adt_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD69', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cd279_adt_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD279', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cd4_adt_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD4', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cd8_adt_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD8', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cd45ra_adt_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD45RA', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
pdf('/filepath/fig1/fig1_supp/tfh_feature_umap/tcell_cd25_adt_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'adt-CD25', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# E) UMAP embedding of L2 Tfh-like clusters vs grouped nonTfh clusters ----

# Grouping nonTfh clusters
tfh_clust_vs_nontfh_group <- ifelse(l2_tcell_clust_names %in% c('Tfh GC','Tfh IL10','CD4 Tcm/fh'), l2_tcell_clust_names, "nonTfh")
names(tfh_clust_vs_nontfh_group) <- names(l2_tcell_clust_names)
l2_teaseq_tcell_obj$tfh_clust_vs_nontfh_group <- l2_teaseq_tcell_obj$wnn.tcell.subcluster.harmony
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
l2_teaseq_tcell_obj <- RenameIdents(l2_teaseq_tcell_obj, tfh_clust_vs_nontfh_group)
l2_teaseq_tcell_obj$tfh_clust_vs_nontfh_group <- Idents(l2_teaseq_tcell_obj)
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
table(l2_teaseq_tcell_obj$tfh_clust_vs_nontfh_group)

# L2 3WNN UMAP without cluster labels - added manually in Adobe Illustrator
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
l2_tfh_vs_nontfh_dimplot <- DimPlot(l2_teaseq_tcell_obj, reduction = "umap.wnn.harmony", pt.size = 1.5, 
                                    group.by = 'tfh_clust_vs_nontfh_group', order = c("CD4 Tcm/fh","Tfh IL10","Tfh GC", 'nonTfh'),
                                    cols = tfh_clust_vs_nontfh_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/l2_tcell_wnn_umap_clust_nolabs.pdf', width = 8, height = 8)
l2_tfh_vs_nontfh_dimplot
dev.off()

# L2 3WNN UMAP with cluster labels for reference
tcell_wnn_umap_coords <- FetchData(l2_teaseq_tcell_obj, vars = c("harmonywnnUMAP_1", "harmonywnnUMAP_2")) %>%
  mutate(tfh_clust_vs_nontfh_group = Idents(l2_teaseq_tcell_obj),
         fill_col = tfh_clust_vs_nontfh_cols[as.character(tfh_clust_vs_nontfh_group)])
l2_wnn_umap_lab_df <- tcell_wnn_umap_coords %>% 
  group_by(tfh_clust_vs_nontfh_group, fill_col) %>% 
  summarize(UMAP_1 = median(harmonywnnUMAP_1), UMAP_2 = median(harmonywnnUMAP_2), .groups = 'drop')
l2_tfh_vs_nontfh_dimplot_lab <- l2_tfh_vs_nontfh_dimplot +
  new_scale_color() + 
  geom_label_repel(
    data = l2_wnn_umap_lab_df,
    aes(x = UMAP_1, y = UMAP_2, label = tfh_clust_vs_nontfh_group, fill = fill_col, color = 'black'),
    size = 6,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines"),
    force = 1,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  NoLegend()
pdf('/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/l2_tfh_vs_nontfh_dimplot_lab.pdf', width = 8, height = 8)
l2_tfh_vs_nontfh_dimplot_lab
dev.off()

# Finding markers of each L2 Tfh cluster and grouped nonTfh clusters across assays ----

# Join layers to compute differentially expressed features
l2_teaseq_tcell_obj <- JoinLayers(l2_teaseq_tcell_obj, assay = 'RNA')
l2_teaseq_tcell_obj <- JoinLayers(l2_teaseq_tcell_obj, assay = 'ADT')

# Find differentially expressed features across assays - ATAC peaks, ATAC GeneActivity, ATAC chromVAR, RNA, RNA SCENIC, and ADT
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
l2_tfh_vs_nontfh_atac_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ATAC', logfc.threshold = 0, min.pct = 0)
l2_tfh_vs_nontfh_act_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ACT',  logfc.threshold = 0, min.pct = 0)
l2_tfh_vs_nontfh_rna_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'RNA', logfc.threshold = 0, min.pct = 0)
l2_tfh_vs_nontfh_adt_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ADT',  logfc.threshold = 0, min.pct = 0)
l2_tfh_vs_nontfh_chromvar_markers <- FindAllMarkers(
  object = l2_teaseq_tcell_obj,
  mean.fxn = rowMeans, # use row means function to find cluster differences as chromVAR Z-scores are not log-normalized like Seurat RNA data
  fc.name = "mean_diff_z_score", # using default parameters
  assay = 'chromvar',
  logfc.threshold = 0, min.pct = 0
)
l2_tfh_vs_nontfh_scenic_markers <- FindAllMarkers(
  object = l2_teaseq_tcell_obj,
  assay = 'SCENIC',
  mean.fxn = rowMeans, # AUC scores, not log-transformed count data
  fc.name = "mean_diff_AUC",
  min.pct = 0, logfc.threshold = 0 # adjust default filters as the mean difference in AUC may be lower than 0.1 
)

# Rename column names for each output FindAllMarkers dataframe for clarity and export ----

# TEAseq 3WNN L3 Tfh-like Cluster Markers - ATAC (peaks)
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
dap_names <- l2_tfh_vs_nontfh_atac_markers$gene
dap_gr <- StringToGRanges(dap_names)
dap_closest_gene <- ClosestFeature(
  object = l2_teaseq_tcell_obj,
  regions    = dap_gr,
  annotation = Annotation(l2_teaseq_tcell_obj)
)
l2_tfh_vs_nontfh_atac_markers <- l2_tfh_vs_nontfh_atac_markers %>%
  mutate(
    closest_gene_symbol = dap_closest_gene$gene_name,
    closest_ensembl_id = dap_closest_gene$gene_id
  )
l2_tfh_vs_nontfh_atac_markers_export <- l2_tfh_vs_nontfh_atac_markers %>%
  rename(
    l2_tfh_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    atac_peak = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_tfh_vs_nontfh_atac_markers_export <- l2_tfh_vs_nontfh_atac_markers_export %>% arrange(l2_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tfh_vs_nontfh_atac_markers_export <- l2_tfh_vs_nontfh_atac_markers_export %>% relocate(c('l2_tfh_cluster','atac_peak','closest_gene_symbol','closest_ensembl_id','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_atac_markers_export, '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/l2_tfh_vs_nontfh_atac_peak_markers_export.rds')

# L2 ATAC GeneActivity spreadsheet export
l2_tfh_vs_nontfh_atac_geneactivity_markers_export <- l2_tfh_vs_nontfh_act_markers %>%
  rename(
    avg_log2fc_clust_vs_rest = avg_log2FC,
    l2_tfh_cluster = cluster,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_nonzero_auc_others = pct.2
  )
l2_tfh_vs_nontfh_atac_geneactivity_markers_export <- l2_tfh_vs_nontfh_atac_geneactivity_markers_export %>% arrange(l2_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tfh_vs_nontfh_atac_geneactivity_markers_export <- l2_tfh_vs_nontfh_atac_geneactivity_markers_export %>% relocate(c('l2_tfh_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_atac_geneactivity_markers_export, '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/l2_tfh_vs_nontfh_atac_geneactivity_markers_export.rds')

# TEAseq 3WNN Each L2 Tfh Cluster vs nonTfh Group Markers - chromVAR (ATAC-based TF motif deviation Z-scores)
pfm <- getMatrixSet( # function to get common TF name for each position weight matrix ID
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
id_to_tf <- sapply(pfm, function(x) name(x))
names(id_to_tf) <- sapply(pfm, ID)
l2_tfh_vs_nontfh_chromvar_markers$tf_name <- id_to_tf[l2_tfh_vs_nontfh_chromvar_markers$gene]
l2_tfh_vs_nontfh_chromvar_markers_export <- l2_tfh_vs_nontfh_chromvar_markers %>%
  rename(
    l2_tfh_cluster = cluster,
    jaspar_pfm_id = gene,
    p_val_raw = p_val,
    pct_nonzero_deviation_in_clust = pct.1, # rename, as nonzero chromVAR Z-scores do not mean 'positive' as with RNA expression 
    pct_nonzero_deviation_in_others = pct.2
  )
l2_tfh_vs_nontfh_chromvar_markers_export <- l2_tfh_vs_nontfh_chromvar_markers_export %>% arrange(l2_tfh_cluster, desc(mean_diff_z_score), p_val_adj)
l2_tfh_vs_nontfh_chromvar_markers_export <- l2_tfh_vs_nontfh_chromvar_markers_export %>% relocate(c('l2_tfh_cluster','tf_name','jaspar_pfm_id','mean_diff_z_score','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_chromvar_markers_export, '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/l2_tfh_vs_nontfh_chromvar_markers_export.rds')

# TEAseq 3WNN Each L2 Tfh Cluster vs nonTfh Group Markers - RNA (log-normalized counts)
l2_tfh_vs_nontfh_rna_markers_export <- l2_tfh_vs_nontfh_rna_markers %>%
  rename(
    avg_log2fc_clust_vs_rest = avg_log2FC,
    l2_tfh_cluster = cluster,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_nonzero_auc_others = pct.2
  )
l2_tfh_vs_nontfh_rna_markers_export <- l2_tfh_vs_nontfh_rna_markers_export %>% arrange(l2_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tfh_vs_nontfh_rna_markers_export <- l2_tfh_vs_nontfh_rna_markers_export %>% relocate(c('l2_tfh_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_rna_markers_export, '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/l2_tfh_vs_nontfh_rna_markers_export.rds')

# TEAseq 3WNN Each L2 Tfh Cluster vs nonTfh Group Markers - SCENIC (RNA-based TF regulon AUC scores)
l2_tfh_vs_nontfh_scenic_markers_export <- l2_tfh_vs_nontfh_scenic_markers %>%
  rename(
    l2_tfh_cluster = cluster,
    scenic_regulon = gene,
    p_val_raw = p_val,
    pct_nonzero_auc_in_clust = pct.1,
    pct_nonzero_auc_others = pct.2
  )
l2_tfh_vs_nontfh_scenic_markers_export <- l2_tfh_vs_nontfh_scenic_markers_export %>% arrange(l2_tfh_cluster, desc(mean_diff_AUC), p_val_adj)
l2_tfh_vs_nontfh_scenic_markers_export <- l2_tfh_vs_nontfh_scenic_markers_export %>% relocate(c('l2_tfh_cluster','scenic_regulon','mean_diff_AUC','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_scenic_markers_export, '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/l2_tfh_vs_nontfh_scenic_markers_export.rds')

# TEAseq 3WNN L3 Tfh-like Cluster Markers - ADT (CLR-normalized count data)
l2_tfh_vs_nontfh_adt_markers_export <- l2_tfh_vs_nontfh_adt_markers %>%
  rename(
    l2_tfh_cluster = cluster,
    adt_epitope = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_tfh_vs_nontfh_adt_markers_export <- l2_tfh_vs_nontfh_adt_markers_export %>% arrange(l2_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tfh_vs_nontfh_adt_markers_export$adt_epitope <- gsub("^adt-", "", l2_tfh_vs_nontfh_adt_markers_export$adt_epitope)
l2_tfh_vs_nontfh_adt_markers_export <- l2_tfh_vs_nontfh_adt_markers_export %>% relocate(c('l2_tfh_cluster','adt_epitope','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_adt_markers_export, '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/l2_tfh_vs_nontfh_adt_markers_export.rds')

# Import marker .rds files and export as separate sheets of single XLSX spreadsheet
clust_marker_df_dir   <- '/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df # prep_clust_markers(df)
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/fig1/fig1_supp/l2_tcell_subclust_annot/l2_teaseq_each_tfh_clust_vs_nontfh_multimodal_markers/teaseq_l2_3wnn_each_tfh_clust_vs_nontfh_markers_all_assays.xlsx")  

# E) Identifying core ATAC features enriched across L2 Tfh-like clusters versus grouped nonTfh clusters ----

# Find ATAC DAPs
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
l2_tfh_vs_nontfh_atac_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ATAC', logfc.threshold = 0, min.pct = 0)

# FIltering for peaks enriched in each group with Bonferroni-adjusted Wilcoxon P < 0.05
ctfh_daps <- l2_tfh_vs_nontfh_atac_markers %>% filter(cluster == 'CD4 Tcm/fh') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) 
il10tfh_daps <- l2_tfh_vs_nontfh_atac_markers %>% filter(cluster == 'Tfh IL10') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) 
gctfh_daps <- l2_tfh_vs_nontfh_atac_markers %>% filter(cluster == 'Tfh GC') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05)
nontfh_daps <-  l2_tfh_vs_nontfh_atac_markers %>% filter(cluster == 'nonTfh') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05)

# Create list of DAP sets from each L2 Tfh and nonTfh subcluster
tfh_compare_dap_sets <- list(
  'CD4 Tcm/fh' = ctfh_daps$gene,
  'Tfh IL10' = il10tfh_daps$gene,
  'Tfh GC' = gctfh_daps$gene,
  'nonTfh' = nontfh_daps$gene
)

# Get cluster DAP set names
tfh_compare_dap_set_names <- names(tfh_compare_dap_sets)

# Get all possible groupings of L2 Tfh and nonTfh
tfh_compare_combo_list <- unlist(lapply(1:length(tfh_compare_dap_set_names),
                                        function(k) combn(tfh_compare_dap_set_names, k, simplify = FALSE)),
                                 recursive = FALSE)

# For each grouping of clusters, find shared DAP distinct from rest of clusters outside grouping
tfh_compare_combo_daps <- lapply(tfh_compare_combo_list, function(combo){
  in_all <- Reduce(intersect, tfh_compare_dap_sets[combo]) # Find DAP common to all clusters of given grouping
  others <- setdiff(tfh_compare_dap_set_names, combo) # find clusters not in given grouping
  if(length(others)) in_all <- setdiff(in_all, unlist(tfh_compare_dap_sets[others])) # Remove DAP found in clusters outside given grouping
  in_all
})
names(tfh_compare_combo_daps) <- sapply(tfh_compare_combo_list, paste, collapse = "&") # Name set comparisons

# Get DAP found only in Tfh clusters versus nonTfh
tfh_daps_df <- as.data.frame(tfh_compare_combo_daps$`CD4 Tcm/fh&Tfh IL10&Tfh GC`) 

# Find closest annotated gene to each DAP
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
tfh_daps_closest_feat <- ClosestFeature( 
  object = l2_teaseq_tcell_obj,
  regions    = tfh_daps_df$`tfh_compare_combo_daps$\`CD4 Tcm/fh&Tfh IL10&Tfh GC\``,
  annotation = Annotation(l2_teaseq_tcell_obj)
)
tfh_daps_closest_feat <- tfh_daps_closest_feat %>% rename(tfh_dap = query_region) %>% select(tfh_dap, everything())
rownames(tfh_daps_closest_feat) <- tfh_daps_closest_feat$tfh_dap
View(tfh_daps_closest_feat)

# Save list of core Tfh like DAP - 415 total unique DAP, near 319 unique gene loci
write_xlsx(tfh_daps_closest_feat, '/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/tfh_specific_peaks_annotated.xlsx')

# E) Identifying core RNA features enriched across L2 Tfh-like clusters versus grouped nonTfh clusters ----

# Find DEGs between each L2 Tfh cluster versus grouped nonTfh clusters
ctfh_degs <- l2_tfh_vs_nontfh_rna_markers %>% filter(cluster == 'CD4 Tcm/fh') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) 
il10tfh_degs <- l2_tfh_vs_nontfh_rna_markers %>% filter(cluster == 'Tfh IL10') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) 
gctfh_degs <- l2_tfh_vs_nontfh_rna_markers %>% filter(cluster == 'Tfh GC') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05)
nontfh_degs <-  l2_tfh_vs_nontfh_rna_markers %>% filter(cluster == 'nonTfh') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05)

# Create DEG lists for Tfh-like and nonTfh subclusters
tfh_compare_deg_sets <- list(
  'CD4 Tcm/fh' = ctfh_degs$gene,
  'Tfh IL10' = il10tfh_degs$gene,
  'Tfh GC' = gctfh_degs$gene,
  'nonTfh' = nontfh_degs$gene
)

# Get DEG list cluster names
tfh_compare_deg_set_names <- names(tfh_compare_deg_sets)

# Get all possible groupings of clusters
tfh_deg_compare_combo_list <- unlist(lapply(1:length(tfh_compare_deg_set_names),
                                            function(k) combn(tfh_compare_deg_set_names, k, simplify = FALSE)),
                                     recursive = FALSE)

# For each combination of clusters, find shared DEG distinct from rest of clusters outside group
tfh_compare_combo_degs <- lapply(tfh_deg_compare_combo_list, function(combo){
  in_all <- Reduce(intersect, tfh_compare_deg_sets[combo]) # get DEG common to all clusters in grouping
  others <- setdiff(tfh_compare_deg_set_names, combo) # find clusters not in selected grouping
  if(length(others)) in_all <- setdiff(in_all, unlist(tfh_compare_deg_sets[others])) # remove DEG that appear in any other cluster grouping
  in_all # return resulting DEG list
})
names(tfh_compare_combo_degs) <- sapply(tfh_deg_compare_combo_list, paste, collapse = "&") # Name comparisons

# Get DEG found only in Tfh-like clusters
tfh_degs_df <- as.data.frame(tfh_compare_combo_degs$`CD4 Tcm/fh&Tfh IL10&Tfh GC`)
colnames(tfh_degs_df) <- 'tfh_deg'
rownames(tfh_degs_df) <- tfh_degs_df$tfh_deg

# Save DEG list
write_xlsx(tfh_degs_df, '/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/tfh_specific_degs.xlsx')

# E) Finding core multimodal ATAC and RNA features enriched across L2 Tfh-like clusters versus grouped nonTfh clusters ----

# Get unique gene names in core ATAC and RNA feature sets
daps_genes <- unique(tfh_daps_closest_feat$gene_name)
length(daps_genes) # 319 genes, some redundancy in nearest genes of the 415 DAPs
degs_genes <- unique(tfh_degs_df$tfh_deg)
length(degs_genes) # 225 genes, all unique as expected

# Find genes present in both DAP and DEG sets
tfh_core_feat <- as.data.frame(intersect(degs_genes, daps_genes))
dim(tfh_core_feat) # 64 multimodal features

# Save list of core L2 Tfh group multimodal features
write_xlsx(tfh_core_feat, '/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/tfh_core_feat.xlsx')

# Find DEG features that were not DAP
genes_in_degs_not_daps <- as.data.frame(setdiff(degs_genes, daps_genes))
dim(genes_in_degs_not_daps) # 161 total
write_xlsx(genes_in_degs_not_daps, '/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/core_tfhgenes_in_degs_not_daps.xlsx')

# Find DAP features that were not DEG
genes_in_daps_not_degs <- as.data.frame(setdiff(daps_genes, degs_genes))
dim(genes_in_daps_not_degs) # 255 total
write_xlsx(genes_in_daps_not_degs, '/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/core_tfh_genes_in_daps_not_degs.xlsx')

# E) Venn diagram visualization of overlaping core ATAC versus RNA features of L2 Tfh-like cells

# Create Venn diagram without labels - added manually in Illustrator
venn_list <- setNames(list(daps_genes, degs_genes), c("", ""))
tfh_core_feat_venn <- ggVennDiagram(
  venn_list,
  label_alpha = 0,
  label_size  = 0
  ) +
  scale_fill_gradient(low = "white", high = "#e17794") +
  theme(legend.position = "none")

# Save plot
pdf('/filepath/fig1/fig1_supp/l2_tcell_tfh_core_dap_deg/tfh_core_feat_venn_nolabs.pdf', width = 6, height = 5)
tfh_core_feat_venn
dev.off()