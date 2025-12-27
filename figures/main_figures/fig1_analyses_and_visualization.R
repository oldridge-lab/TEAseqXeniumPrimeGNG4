# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. 1
# Fig. 1 - Trimodal analysis resolves distinct Tfh-like states among tonsil and peripheral blood mononuclear cells.

# Set up working environment ----

# Set working directory to Fig 1
setwd('/filepath/fig1/')

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
library(ggrepel) # 0.9.5

# Import Seurat object from TEAseq Data Preprocessing Step 13 (trimodal dimensionality reduction and L1 3WNN clustering of all tonsil and peripheral blood mononuclear cells, including Harmony integration across donors)
l1_teaseq_obj <- readRDS('/filepath/step13_bulk_harmony/l1_teaseq_obj.rds')

# Import Seurat object from TEAseq Data Preprocessing Step 14 (trimodal dimensionality reduction and L2 3WNN subclustering of T cells from L1 object, including Harmony integration across donors)
l2_teaseq_tcell_obj <- readRDS('/filepath/step14_tcell_subcluster/l2_teaseq_tcell_obj.rds')

# Specify colors for Level 1 3WNN clusters
l1_clust_cols <- c(
  "CD4 Tn" = "#8ab1cc",       
  "NBC" = "#657ab0",        
  "CD4 Tm" = "#94a890",      
  "CD8 Tm/Inn" = "#add5df",
  "Tfh-like" = "#e49a78",         
  "MBC" = "#d4b5e4",         
  "CD8 Tn" = "#c49980",      
  "NK" = "#728782",      
  "Treg" = "#c9a3c1",        
  "GCB" = "#d66d97",         
  "Myl" = "#b6b7ba",         
  "CD4 Tn Ribo" = "#e0c1b4",        
  "Myl CD16" = "#e1d4c7",    
  "CD8 UTC" = "#d699ac",  
  "ASC" = "#dfcc78"          
)

# Specify colors for each tissue
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# Fig 1A - TEAseq experimental design schematic ----

# Designed in BioRender and exported as PDF

# Fig 1B - 3WNN L1 clustering and UMAP visualization ----

# Set 3WNN L1 cluster annotations
Idents(l1_teaseq_obj) <- "wnn.bulk.cluster.harmony" # From 0.2 resolution L1 3WNN clustering, refer to Step 13 code
l1_wnn_clust_names <- c(
  "0" = "CD4 Tn",
  "1" = "NBC",
  "2" = "CD4 Tm",
  "3" = "CD8 Tm/Inn",
  "4" = "Tfh-like",
  "5" = "MBC",
  "6" = "CD8 Tn",
  "7" = "NK",
  "8" = "Treg",
  "9" = "GCB",
  "10" = "Myl",
  "11" = "CD4 Tn Ribo",
  "12" = "Myl CD16",
  "13" = "CD8 UTC",
  "14" = "ASC"
)

# Rename 15 coarse 3WNN cluster numbers to assigned names
l1_teaseq_obj$l1_wnn_annot <- l1_teaseq_obj$wnn.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_wnn_clust_names)
l1_teaseq_obj$l1_wnn_annot <- Idents(l1_teaseq_obj)

# L1 3WNN UMAP including cluster labels
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
l1_wnn_umap_clust_dimplot <- DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", pt.size = 0.5, group.by = 'l1_wnn_annot', cols = l1_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
l1_wnn_umap_coords <- FetchData(l1_teaseq_obj, vars = c("harmonywnnUMAP_1", "harmonywnnUMAP_2")) %>%
  mutate(l1_wnn_annot = Idents(l1_teaseq_obj),
         fill_col = l1_clust_cols[as.character(l1_wnn_annot)])
l1_wnn_umap_lab_df <- l1_wnn_umap_coords %>% 
  group_by(l1_wnn_annot, fill_col) %>% 
  summarize(UMAP_1 = median(harmonywnnUMAP_1), UMAP_2 = median(harmonywnnUMAP_2), .groups = 'drop')
l1_wnn_umap_clust_lab <- l1_wnn_umap_clust_dimplot +
  new_scale_color() + 
  geom_label_repel(
    data = l1_wnn_umap_lab_df,
    aes(x = UMAP_1, y = UMAP_2, label = l1_wnn_annot, fill = fill_col, color = 'black'),
    size = 6,
    box.padding = unit(0.7, "lines"),
    point.padding = unit(1.5, "lines"),
    label.padding = unit(0.35, "lines"),
    force = 2,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  NoLegend()
pdf("/filepath/fig1b/l1_wnn_umap_clust_lab.pdf", width = 6, height = 6)
l1_wnn_umap_clust_lab
dev.off()

# L1 3WNN UMAP without cluster labels - applied manually in Adobe Illustrator
pdf("/filepath/fig1b/l1_wnn_umap_clust_no_lab.pdf", width = 6, height = 6)
l1_wnn_umap_clust_dimplot
dev.off()

# Fig 1C - L1 3WNN UMAP split by PBMC vs Tonsil sample origin ----

l1_wnn_umap_tissue <- DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = FALSE, pt.size = 0.5, repel = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf("/filepath/fig1c/l1_wnn_umap_tissue.pdf", width = 6, height = 6)
l1_wnn_umap_tissue
dev.off()

# Fig 1D - Stacked barplot of donor proportion per L1 3WNN cluster, normalized to donor cell total and rescaled to sum to 1 ----

# Determine normalized and rescaled proportions per cluster only using unenriched samples (i.e. not CD4-enriched samples)
Idents(l1_teaseq_obj) <- 'hto.sort'
l1_unenr_samp_md <- subset(l1_teaseq_obj, idents = c('P1-Bulk','P2-Bulk','P3-Bulk','P4-Bulk','T1-Bulk','T2-Bulk','T3-Bulk','T4-Bulk'))@meta.data 
donor_totals <- l1_unenr_samp_md %>%
  group_by(hto.donor) %>%
  summarise(total = n(), .groups = "drop") # number of cells from each donor in each cluster
norm_prop_df <- l1_unenr_samp_md %>%
  group_by(l1_wnn_annot, hto.donor) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(donor_totals, by = "hto.donor") %>%
  mutate(raw_fraction = count / total) # within donor proportion of cells corresponding to given cluster
norm_prop_df <- norm_prop_df %>%
  group_by(l1_wnn_annot) %>%
  mutate(sum_raw = sum(raw_fraction),
         normalized_fraction = raw_fraction / sum_raw) %>% # rescale proportions to sum to 1
  ungroup()

# Find clusters enriched in tonsil samples, and sort clusters from low to high for barplot visualization
t_proportions <- norm_prop_df %>%
  filter(hto.donor %in% paste0("T", 1:4)) %>%
  group_by(l1_wnn_annot) %>%
  summarise(t_sum = sum(normalized_fraction)) %>%
  complete(l1_wnn_annot = levels(norm_prop_df$l1_wnn_annot), fill = list(t_sum = 0)) %>%
  arrange(t_sum)
barplot_cluster_order <- t_proportions$l1_wnn_annot
norm_prop_df <- norm_prop_df %>% mutate(l1_wnn_annot = factor(l1_wnn_annot, levels = barplot_cluster_order))

# Specify colors for each tissue and donor
pbmc_colors <- colorRampPalette(c("steelblue4", "steelblue1"))(4)
tonsil_colors <- colorRampPalette(c("darkorange3", "orange"))(4)
names(pbmc_colors) <- paste0("P", 1:4)
names(tonsil_colors) <- paste0("T", 1:4)
donor_colors <- c(pbmc_colors, tonsil_colors)

# Assemble barplot
l1_clust_prop_norm_barplot <- ggplot(norm_prop_df, aes(x = factor(l1_wnn_annot), y = normalized_fraction, fill = hto.donor)) +
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
l1_clust_prop_norm_barplot
ggsave("/filepath/fig1d/l1_clust_prop_norm_barplot.pdf", plot = l1_clust_prop_norm_barplot, width = 12, height = 8)

# Fig 1E - Ternary plot of WNN modality weights for centroid of each L1 cluster ----
l1_teaseq_obj.md <- l1_teaseq_obj@meta.data
l1_teaseq_obj.md <- l1_teaseq_obj.md %>%
  rename(RNA = RNA.weight, ATAC = ATAC.weight, ADT = ADT.weight)
bulk.tern.centroids <- l1_teaseq_obj.md %>%
  group_by(l1_wnn_annot) %>%
  summarise(
    ATAC = mean(ATAC, na.rm = TRUE),
    RNA = mean(RNA, na.rm = TRUE),
    ADT = mean(ADT, na.rm = TRUE)
    
  )
tern_point_order <- c('NBC','MBC','GCB','ASC','Myl','CD8 Tn','CD8 Tm/Inn','CD8 UTC','NK','Myl CD16','CD4 Tn','CD4 Tn Ribo','CD4 Tm','Treg','Tfh-like')
bulk.tern.centroids.sorted <- bulk.tern.centroids %>%
  mutate(l1_wnn_annot = factor(l1_wnn_annot, levels = tern_point_order)) %>%
  arrange(l1_wnn_annot)
options(tern.arrow = arrow(type = "open", length = unit(0.5, "cm")))
l1_clust_mod_weights_tern_legend <- ggtern(data = bulk.tern.centroids.sorted, aes(x = ATAC, y = ADT, z = RNA)) +
  geom_point(aes(fill = l1_wnn_annot), size = 8, shape = 21, color = 'black', stroke = 1) +
  scale_fill_manual(values = l1_clust_cols, name = 'l1_wnn_annot') +
  labs(color = "Cluster") +
  theme_bw(base_family = 'sans') +
  theme_hidegrid_minor() +
  theme_showarrows() +
  theme_arrowlarge() +
  geom_confidence_tern(color = 'black') +
  theme_hidetitles() +
  theme(
    text               = element_text(family = "sans"),
    axis.title         = element_text(family = "sans", color = 'black'),
    axis.text          = element_text(family = "sans", color = 'black'),
    legend.text        = element_text(size = 10, family = "sans", color = 'black'),
    legend.title       = element_blank(),
    tern.axis.text = element_text(family = 'sans', size= 20),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "black"),
    panel.grid.minor = element_line(color = "black"),
    tern.axis.arrow = element_line(linewidth = 3, color = 'black'),
    tern.axis.arrow.sep  = 0.15,
    tern.axis.vshift     = 0.05,
    legend.position = 'bottom'
  ) +
  guides(fill = guide_legend(ncol = 5, byrow = TRUE))
l1_clust_mod_weights_tern_legend
ggsave("/filepath/fig1e/l1_clust_mod_weights_tern_legend.pdf", plot = l1_clust_mod_weights_tern_legend, width = 8, height = 6)

# Fig 1F - 3WNN L1 cluster annotations applied to unimodal ATAC vs RNA vs ADT UMAP embeddings----

# ATAC UMAP L1
l1_atac_umap_clust_lab <- DimPlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = FALSE, pt.size = 0.5, group.by = 'l1_wnn_annot', cols = l1_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1f/l1_atac_umap_clust_lab.pdf', width = 6, height = 6)
l1_atac_umap_clust_lab
dev.off()

# RNA UMAP L1
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
l1_rna_umap_clust_lab_rot <- DimPlot(l1_teaseq_obj, reduction = "umap.rna.harmony_rot", label = FALSE, pt.size = 0.5, group.by = 'l1_wnn_annot', cols = l1_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1f/l1_rna_umap_clust_lab_rot.pdf', width = 6, height = 6)
l1_rna_umap_clust_lab_rot
dev.off()

# ADT UMAP L1
l1_adt_umap_clust_lab <- DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = FALSE, pt.size = 0.5, group.by = 'l1_wnn_annot', cols = l1_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1f/l1_adt_umap_clust_lab.pdf', width = 6, height = 6)
l1_adt_umap_clust_lab
dev.off()

# Fig 1G - Barplot of adjusted mutual information of 15x L1 cluster identities between unimodal analyses ----

# Use scikit-learn Python package to determine clustering similarity metrics
library(reticulate) # 1.37.0
py_install("scikit-learn", pip = TRUE)
sklearn <- import("sklearn.metrics")
py_version() # 3.9
pkg_res <- import("pkg_resources", convert = FALSE)
pkg_res$get_distribution("scikit-learn")$version # 1.6.1

# Get numerical cluster assignments (refer to 'figS5_L1_unimodal_ATAC_RNA_ADT_clust_annotation.R' for full clustering information)
cluster_numb_RNA <- l1_teaseq_obj$rna.bulk.cluster.harmony
cluster_numb_ADT <- l1_teaseq_obj$adt.bulk.cluster.harmony
cluster_numb_ATAC <- l1_teaseq_obj$atac.bulk.cluster.harmony

# Compute AMI between unimodal clustering sets
allclust_ami_RNA_ADT <- sklearn$adjusted_mutual_info_score(cluster_numb_RNA, cluster_numb_ADT)
allclust_ami_RNA_ATAC <- sklearn$adjusted_mutual_info_score(cluster_numb_RNA, cluster_numb_ATAC)
allclust_ami_ADT_ATAC <- sklearn$adjusted_mutual_info_score(cluster_numb_ADT, cluster_numb_ATAC)

# Assemble dataframe
l1_clust_modality_ami_df <- data.frame(
  comparison = c("ATAC vs RNA","RNA vs ADT","ATAC vs ADT"),
  metric = rep("AMI", times = 3),
  value = c(0.6181634929906691,
            0.7004975674278031,
            0.5779202550584999)
)

# Set plot order
l1_clust_modality_ami_df$comparison <- factor(l1_clust_modality_ami_df$comparison, levels = c("ATAC vs RNA","RNA vs ADT", "ATAC vs ADT"))
l1_clust_modality_ami_df <- l1_clust_modality_ami_df %>% mutate(value_rounded = round(value, digits = 2))
modality_comp_cols <- c('ATAC vs RNA' = 'grey', 'RNA vs ADT'= 'grey', 'ATAC vs ADT' = 'grey')

# Create barplot
l1_clust_modality_ami_barplot <- ggplot(l1_clust_modality_ami_df, aes(x = comparison, y = value_rounded, fill = comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = 'black', linewidth = 0.75) +
  scale_fill_manual(values = modality_comp_cols) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  labs(title = "L1 Clustering Similarity",
       x = "comparison",
       y = "Adjusted Mutual Information",
       fill = "metric") +
  theme_classic(base_size = 20, 'sans') +
  theme(axis.text.x = element_text(angle = 45, 'sans', hjust = 1, margin = margin(t = 5), size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5),  size = 20, 'sans'),
        plot.title = element_blank()) + 
  geom_text(aes(label = formatC(value_rounded, format = "f", digits = 2)),
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 7) +
  NoLegend()

# Remove x-axis labels - applied in Illustrator
pdf("/filepath/fig1g/l1_clust_modality_ami_barplot_no_lab_grey.pdf", width = 5, height = 4.5)
l1_clust_modality_ami_barplot + theme(
  axis.ticks = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 18, 'sans'),
  axis.title.y = element_text(margin = margin(r = 5), size = 22, 'sans')
)
dev.off()

# Fig 1H - L2 T cell 3WNN subclustering  ----

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

# Specify L2 TEAseq T cell cluster colors
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

# L2 3WNN UMAP without cluster labels - added manually in Adobe Illustrator
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'
tcell_wnn_umap_clust_dimplot <- DimPlot(l2_teaseq_tcell_obj, reduction = "umap.wnn.harmony", pt.size = 1.5, group.by = 'l2_tcell_wnn_annot', cols = tcell_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1h/l2_tcell_wnn_umap_clust_nolabs.pdf', width = 8, height = 8)
tcell_wnn_umap_clust_dimplot
dev.off()

# L2 3WNN UMAP with cluster labels for reference
tcell_wnn_umap_coords <- FetchData(l2_teaseq_tcell_obj, vars = c("harmonywnnUMAP_1", "harmonywnnUMAP_2")) %>%
  mutate(l2_tcell_wnn_annot = Idents(l2_teaseq_tcell_obj),
         fill_col = tcell_clust_cols[as.character(l2_tcell_wnn_annot)])
l1_wnn_umap_lab_df <- tcell_wnn_umap_coords %>% 
  group_by(l2_tcell_wnn_annot, fill_col) %>% 
  summarize(UMAP_1 = median(harmonywnnUMAP_1), UMAP_2 = median(harmonywnnUMAP_2), .groups = 'drop')
tcell_wnn_umap_clust_lab <- tcell_wnn_umap_clust_dimplot +
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
  scale_fill_identity() +
  scale_color_identity() +
  NoLegend()
pdf('/filepath/fig1h/l2_tcell_wnn_umap_clust_lab.pdf', width = 8, height = 8)
tcell_wnn_umap_clust_lab
dev.off()

# L2 3WNN T Cell UMAP split by tissue
tcell_wnn_umap_tissue <- DimPlot(l2_teaseq_tcell_obj, reduction = "umap.wnn.harmony", label = FALSE, label.box = TRUE, label.size = 4, pt.size = 1.5, repel = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig1h/tcell_wnn_umap_tissue.pdf', width = 8, height = 8)
tcell_wnn_umap_tissue
dev.off()

# Visualization of L2 T cell subcluster features on 3WNN UMAP

DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
pdf('/filepath/fig1h/tfh_feature_umap/tcell_cxcr5_prom_umap.pdf', width = 6, height = 6) # CXCR5 promoter-like sequence-containing DAP, as annotated in ENCODE
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'chr11-118883255-118884245', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1h/tfh_feature_umap/tcell_cxcr5_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'CXCR5', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1h/tfh_feature_umap/tcell_tox2_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'TOX2', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
pdf('/filepath/fig1h/tfh_feature_umap/tcell_il21_rna_umap.pdf', width = 6, height = 6)
do_FeaturePlot(l2_teaseq_tcell_obj, reduction = 'umap.wnn.harmony', features = 'IL21', use_viridis = FALSE, sequential.palette = "Oranges", order = TRUE, legend.position = 'right', legend.length = 11, border.size = 1.25, pt.size = 1) + coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig 1I - DEP Volcano Plot in Tfh vs nonTfh ----

tfh_clust_vs_nontfh_group <- ifelse(l2_tcell_clust_names %in% c('Tfh GC','Tfh IL10','CD4 Tcm/fh'), l2_tcell_clust_names, "nonTfh")
names(tfh_clust_vs_nontfh_group) <- names(l2_tcell_clust_names)
l2_teaseq_tcell_obj$tfh_clust_vs_nontfh_group <- l2_teaseq_tcell_obj$wnn.tcell.subcluster.harmony
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
l2_teaseq_tcell_obj <- RenameIdents(l2_teaseq_tcell_obj, tfh_clust_vs_nontfh_group)
l2_teaseq_tcell_obj$tfh_clust_vs_nontfh_group <- Idents(l2_teaseq_tcell_obj)
Idents(l2_teaseq_tcell_obj) <- 'tfh_clust_vs_nontfh_group'
table(l2_teaseq_tcell_obj$tfh_clust_vs_nontfh_group)

l2_teaseq_tcell_obj_joined <- JoinLayers(l2_teaseq_tcell_obj, assay = 'ADT')
Idents(l2_teaseq_tcell_obj_joined) <- 'tfh_clust_vs_nontfh_group'
tfh_vs_nontfh_dep <- FindMarkers(l2_teaseq_tcell_obj_joined, ident.1 = c('CD4 Tcm/fh','Tfh IL10', 'Tfh GC'), ident.2 = 'nonTfh', assay = 'ADT',
                                 min.pct = 0, logfc.threshold = 0)
write.csv(tfh_vs_nontfh_dep, '/filepath/fig1i/tfh_vs_nontfh_dep.csv', row.names = FALSE)
saveRDS(tfh_vs_nontfh_dep, '/filepath/fig1i/tfh_vs_nontfh_dep.rds')

# Assemble volcano plot
tfh_vs_nontfh_dep_volc_cols <- ifelse(
  tfh_vs_nontfh_dep$avg_log2FC < -0.3 & tfh_vs_nontfh_dep$p_val < 1e-30,  "#1f77b4",  
  ifelse(
    tfh_vs_nontfh_dep$avg_log2FC >  0.3 & tfh_vs_nontfh_dep$p_val < 1e-30,  "darkorange3",
    "lightgrey"                            
  )
)
# Name the keyvals so they appear in the legend
names(tfh_vs_nontfh_dep_volc_cols)[tfh_vs_nontfh_dep_volc_cols == "lightgrey"] <- "Nonsignificant"
names(tfh_vs_nontfh_dep_volc_cols)[tfh_vs_nontfh_dep_volc_cols == "#1f77b4"]    <- "Up in nonTfh-like Group"
names(tfh_vs_nontfh_dep_volc_cols)[tfh_vs_nontfh_dep_volc_cols == "darkorange3"]    <- "Up in Tfh-like Group"
tfh_vs_nontfh_dep_names_trim <- gsub("^adt-", "", rownames(tfh_vs_nontfh_dep))
rownames(tfh_vs_nontfh_dep) <- tfh_vs_nontfh_dep_names_trim
tfh_vs_nontfh_dep_name_map <- c(
  "CD278" = "ICOS",
  "CD279" = "PD-1",
  "CD185" = "CXCR5",
  "CD304" = "NRP-1",
  "CD11a" = "LFA-1a",
  "CD152" = "CTLA-4",
  "CD305" = "LAIR-1",
  "CD314" = "NKG2D",
  "CD29"  = "ITGB1",
  "CD162" = "PSGL-1",
  "CD58"  = "LFA-3",
  "CD18"  = "LFA-1b",
  "HLA.DR"= "HLA-DR",
  "HLA.DR.DP.DQ" = "MHC-II"
)
tfh_vs_nontfh_dep$adt <- rownames(tfh_vs_nontfh_dep)
tfh_vs_nontfh_dep_rownames <- rownames(tfh_vs_nontfh_dep)
dep_name_index <- which(tfh_vs_nontfh_dep_rownames %in% names(tfh_vs_nontfh_dep_name_map))
tfh_vs_nontfh_dep_rownames[dep_name_index] <- tfh_vs_nontfh_dep_name_map[ tfh_vs_nontfh_dep_rownames[dep_name_index] ]
rownames(tfh_vs_nontfh_dep) <- tfh_vs_nontfh_dep_rownames

# Save ADT list
l2_tfh_vs_nontfh_adt_export <- tfh_vs_nontfh_dep %>%
  rename(
    avg_log2fc_tfh_vs_nontfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_tfh = pct.1,
    pct_pos_in_nontfh = pct.2
  )
l2_tfh_vs_nontfh_adt_export <- l2_tfh_vs_nontfh_adt_export %>% arrange(desc(avg_log2fc_tfh_vs_nontfh), p_val_adj)
l2_tfh_vs_nontfh_adt_export <- l2_tfh_vs_nontfh_adt_export %>% relocate(c('adt','avg_log2fc_tfh_vs_nontfh','p_val_raw','p_val_adj'))
saveRDS(l2_tfh_vs_nontfh_adt_export, '/filepath/fig1i/l2_tfh_vs_nontfh_adt_export.rds')
write_xlsx(l2_tfh_vs_nontfh_adt_export, '/filepath/fig1i/l2_tfh_vs_nontfh_adt_export.xlsx')

# ADT volcano plot without labels - added manually in Adobe Illustrator to improve visualization
tfh_vs_nontfh_dep_vol_labs <- tfh_vs_nontfh_dep %>% filter(p_val < 1e-30) %>% filter(abs(avg_log2FC) > 0.3) %>% arrange(p_val_adj) %>% slice_head(n = 50) %>% rownames()
pdf('/filepath/fig1j/tfh_vs_nontfh_dep_volcano_nolabs.pdf', height = 9.25, width = 8)
EnhancedVolcano(tfh_vs_nontfh_dep,
                lab = rownames(tfh_vs_nontfh_dep),
                selectLab = c(''),
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Tfh vs nonTfh-like DEP',
                drawConnectors = TRUE,
                pCutoff = 1e-30,
                FCcutoff = 0.3,
                pointSize = 5,
                labSize = 5,
                colCustom= tfh_vs_nontfh_dep_volc_cols,
                legendPosition = 'none',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  xlim(-1.81,3.2)
dev.off()


# Volcano plot with ADT labels for reference
pdf('/filepath/fig1j/tfh_vs_nontfh_dep_volcano_labs.pdf', height = 9.25, width = 8)
EnhancedVolcano(tfh_vs_nontfh_dep,
                lab = rownames(tfh_vs_nontfh_dep),
                selectLab = tfh_vs_nontfh_dep_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Tfh vs nonTfh-like DEP',
                drawConnectors = TRUE,
                pCutoff = 1e-30,
                FCcutoff = 0.3,
                pointSize = 5,
                labSize = 5,
                colCustom= tfh_vs_nontfh_dep_volc_cols,
                legendPosition = 'none',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  xlim(-1.81,3.2)
dev.off()

# Fig 1J - CXCR5 Locus Accessibility in T Cell Subclusters ----

# Prepare GRanges for CXCR5 DAPs of interest between L2 T cell clusters
# chr11-118883255-118884245 = CXCR5 promoter-like sequence-containing DAP enriched in Tfh subclusters
# chr11-118872101-118873015 = CXCR5 distal enhancer-like sequence-containg DAP enriched in Tfh subclusters
tfh_cxcr5_dap_df <- data.frame(
  seqnames = "chr11",
  start = c(118883255, 118872101), # Coordinates for start of CXCR5 PLS- and dELS-containing DAPs, respectively
  end   = c(118884245, 118873015) # Coordinates for end of CXCR5 PLS- and dELS-containing DAPs, respectively
)
tfh_cxcr5_dap_gr <- makeGRangesFromDataFrame(tfh_cxcr5_dap_df)
tfh_cxcr5_dap_gr$color <- c('darkorange3','darkorange3')

# Link CXCR5 peaks to gene expression
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
l2_teaseq_tcell_obj <- RegionStats(l2_teaseq_tcell_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
l2_teaseq_tcell_obj <- JoinLayers(l2_teaseq_tcell_obj, assay = 'RNA')
l2_teaseq_tcell_obj <- LinkPeaks(
  object = l2_teaseq_tcell_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = 'CXCR5')

# Set cluster order for visualization
Idents(l2_teaseq_tcell_obj) <- 'tcell_wnn_annot'
wnn_tcell_dotplot_order <- c("CD4 Tn 1", "CD4 Tn 2", "CD4 Tn 3", "CD4 Tn 4", "CD4 Tn 5", "CD4 Tn 6", "CD4 Tn 7","CD4 Tem", "Th17", "Treg RORgt", "cTreg", "CD4 Tcm/fh", "Tfh IL10", "Tfh GC", "CD4 Inn", "PLZF Inn", "gdT", "CD8 Tn", "CD8 Tm", "CD8 CTL", "CD8 UTC")
Idents(l2_teaseq_tcell_obj) <- factor(
  Idents(l2_teaseq_tcell_obj),
  levels = wnn_tcell_dotplot_order
)

# Assemble CXCR5 coverage plot across clusters
cxcr5_cov_plot <- CoveragePlot(l2_teaseq_tcell_obj, region = 'CXCR5', annotation = TRUE, peaks = TRUE, links = TRUE, region.highlight = tfh_cxcr5_dap_gr, extend.upstream = 12000, extend.downstream = 0, 
                               heights = c(7, 1, 1, 1)) & scale_fill_manual(values = tcell_clust_cols)
cxcr5_cov_plot[[1]][[1]] <- cxcr5_cov_plot[[1]][[1]] + theme(axis.ticks = element_blank()) + labs(y = 'Normalized Peak Signal (0-220)')
cxcr5_cov_plot[[1]][[2]] <- cxcr5_cov_plot[[1]][[2]] + scale_color_manual(values = c('grey','grey'))
cxcr5_cov_plot[[1]][[3]] <- cxcr5_cov_plot[[1]][[3]] + scale_color_manual(values = c('grey','grey','grey','grey','grey','grey','grey','grey'))
cxcr5_cov_plot_links <- cxcr5_cov_plot[[1]][[4]]
cxcr5_cov_plot_links <- cxcr5_cov_plot_links +
  scale_color_gradient(
    low  = "lightgrey",
    high = "darkorange3",
    name = "" 
  ) +
  theme(
    legend.key.size   = unit(0.35, "cm"),
    legend.text       = element_text(size = 10, family = 'sans'),
    legend.title      = element_text(size = 10, family = 'sans')
  )
cxcr5_cov_plot[[1]][[4]] <- cxcr5_cov_plot_links
pdf('/filepath/fig1i/tcell_cxcr5_link_plot.pdf', width = 6, height = 6)
cxcr5_cov_plot
dev.off()

# CXCR5 DotPlot across L2 T cell clusters
DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
Idents(l2_teaseq_tcell_obj) <- factor(
  Idents(l2_teaseq_tcell_obj),
  levels = rev(wnn_tcell_dotplot_order)
)
wnn_tcell_cxcr5_dotplot_feats <- c('CXCR5')
puor_colors <- brewer.pal(11, "PuOr")
wnn_tcell_cxcr5_dotplot <- DotPlot(l2_teaseq_tcell_obj, features = wnn_tcell_cxcr5_dotplot_feats,
                                   cluster.idents = FALSE, cols = 'PuOr', scale = TRUE, dot.scale = 9.75) + 
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",    
        legend.box.just = "center", 
        legend.box = "vertical",             
        legend.title = element_blank(),
        legend.text = element_text(size = 8, family = 'sans'),
        legend.margin = margin(t = 2, b = 2),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 0.1)) +
  guides(
    size = guide_legend(
      title = "",
      order = 1
    ),
    color = guide_colorbar(
      title = "",
      order = 2,
      barwidth = 7,
      barheight = 0.4
    )
  )  + 
  scale_color_gradient2(
    low = puor_colors[11],
    mid = "white",
    high = puor_colors[3],
    midpoint = 0,
    limits = c(-0.5, 2.5)
  )
pdf('/filepath/fig1i/wnn_tcell_cxcr5_dotplot.pdf', height = 8, width = 1.75)
wnn_tcell_cxcr5_dotplot
dev.off()