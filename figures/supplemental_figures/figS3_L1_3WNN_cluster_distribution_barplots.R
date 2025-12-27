# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S3
# Fig. S3 - Distribution of Level 1 cell types in TEAseq dataset across biological and technical variables. (Related to Fig. 1D)

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

# Import Seurat object from TEAseq Data Preprocessing Step 13 (trimodal dimensionality reduction and L1 3WNN clustering of all tonsil and peripheral blood mononuclear cells, including Harmony integration across donors ----
l1_teaseq_obj <- readRDS('/filepath/step13_bulk_harmony/l1_teaseq_obj.rds')

# Order annotation names for plots
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
wnn_dotplot_order <- c('CD4 Tn','CD4 Tn Ribo','CD4 Tm','Treg','Tfh-like','CD8 Tn','CD8 Tm/Inn','CD8 UTC','NK','NBC','MBC','GCB','ASC','Myl','Myl CD16')
Idents(l1_teaseq_obj) <- factor(
  Idents(l1_teaseq_obj),
  levels = wnn_dotplot_order
)

# A) Barplot of L1 3WNN cluster proportions per donor (unenriched samples only) ----

# only including bulk-sorted samples (i.e. excluding CD4-enriched samples, as in Fig 1C)
Idents(l1_teaseq_obj) <- 'hto.sort'
l1_unenr_samp_md <- subset(l1_teaseq_obj, idents = c('P1-Bulk','P2-Bulk','P3-Bulk','P4-Bulk','T1-Bulk','T2-Bulk','T3-Bulk','T4-Bulk'))@meta.data # For calculating donor proportion per cluster, only using samples containing all mononuclear cell types - not the '-CD4' suffix samples enriched for CD4+

# Get cell total per donor
donor_totals <- l1_unenr_samp_md %>%
  group_by(hto.donor) %>%
  summarise(total = n(), .groups = "drop")

# Get proportion of each cluster for each donor
cluster_prop_df <- l1_unenr_samp_md %>%
  group_by(l1_wnn_annot, hto.donor) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(donor_totals, by = "hto.donor") %>%
  mutate(raw_fraction = count / total)

# Make barplot
l1_clust_prop_per_donor_bulk_samples_only_barplot <- ggplot(cluster_prop_df, aes(x = hto.donor, y = raw_fraction, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = l1_clust_cols) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', margin = margin(t = 10)),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_per_donor_bulk_samples_only_barplot
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_per_donor_bulk_samples_only_barplot.pdf", plot = l1_clust_prop_per_donor_bulk_samples_only_barplot, width = 10, height = 6)

# B) Barplot of L1 cluster proportion for each donor, with bulk and sorted samples for each donor combined ----
all_samples_md <- l1_teaseq_obj@meta.data
donor_cell_tot <- all_samples_md %>%
  group_by(hto.donor) %>%
  summarise(total = n(), .groups = "drop")
donor_cell_prop_df <- all_samples_md %>%
  group_by(l1_wnn_annot, hto.donor) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(donor_cell_tot, by = "hto.donor") %>%
  mutate(prop = count / total)
l1_clust_prop_per_donor_bulk_and_enriched_samples_combined_barplot <- ggplot(donor_cell_prop_df, aes(x = hto.donor, y = prop, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = l1_clust_cols) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster",
       x = 'Donor') +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', margin = margin(t = 10)),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_per_donor_bulk_and_enriched_samples_combined_barplot # Unlike Fig 1C, here both bulk-sorted and CD4-enriched samples are used to show cell abundance per donor
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_per_donor_bulk_and_enriched_samples_combined_barplot.pdf", plot = l1_clust_prop_per_donor_bulk_and_enriched_samples_combined_barplot, width = 10, height = 6)

# C) Barplot of L1 cluster proportion for bulk-sorted and CD4-enriched samples separately, rather than grouped per donor ----
all_samples_md <- l1_teaseq_obj@meta.data
sample_cell_tot <- all_samples_md %>%
  group_by(hto.sort) %>%
  summarise(total = n(), .groups = "drop")
sample_cell_prop_df <- all_samples_md %>%
  group_by(l1_wnn_annot, hto.sort) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(sample_cell_tot, by = "hto.sort") %>%
  mutate(prop = count / total)
l1_clust_prop_per_sample_barplot <- ggplot(sample_cell_prop_df, aes(x = hto.sort, y = prop, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = l1_clust_cols) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_per_sample_barplot 
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_per_sample_barplot.pdf", plot = l1_clust_prop_per_sample_barplot, width = 14, height = 7)

# D) Barplot of L1 cluster count per bulk-sorted or CD4-enriched sample separately, rather than grouped per donor ----
l1_clust_count_per_sample_barplot <- ggplot(sample_cell_prop_df, aes(x = hto.sort, y = count, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = l1_clust_cols) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Cells per L1 Cluster",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_count_per_sample_barplot 
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_count_per_sample_barplot.pdf", plot = l1_clust_count_per_sample_barplot, width = 14, height = 7)

# E) Barplot of L1 cluster proportion in CD4-enriched versus bulk cell samples ----
all_samples_md <- l1_teaseq_obj@meta.data
l1_teaseq_obj$cd4_vs_bulk <- ifelse(grepl("-CD4$", l1_teaseq_obj$hto.sort), "CD4-Enriched",
                                    ifelse(grepl("-Bulk$", l1_teaseq_obj$hto.sort), "All Cells", NA))
cd4_vs_bulk_cell_tot <- all_samples_md %>%
  group_by(cd4_vs_bulk) %>%
  summarise(total = n(), .groups = "drop")
cd4_vs_bulk_cell_prop_df <- all_samples_md %>%
  group_by(l1_wnn_annot, cd4_vs_bulk) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(cd4_vs_bulk_cell_tot, by = "cd4_vs_bulk") %>%
  mutate(prop = count / total)
l1_clust_prop_cd4_vs_bulk_barplot <- ggplot(cd4_vs_bulk_cell_prop_df, aes(x = cd4_vs_bulk, y = prop, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = l1_clust_cols) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_cd4_vs_bulk_barplot 
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_cd4_vs_bulk_barplot.pdf", plot = l1_clust_prop_cd4_vs_bulk_barplot, width = 8, height = 8)

# F) Barplot of L1 cluster proportion per tissue, including bulk and CD4-enriched samples together ----
tissue_cell_tot <- all_samples_md %>%
  group_by(hto.tissue) %>%
  summarise(total = n(), .groups = "drop")
tissue_cell_prop_df <- all_samples_md %>%
  group_by(l1_wnn_annot, hto.tissue) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(tissue_cell_tot, by = "hto.tissue") %>%
  mutate(prop = count / total)
l1_clust_prop_per_tissue_barplot <- ggplot(tissue_cell_prop_df, aes(x = hto.tissue, y = prop, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = l1_clust_cols) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_per_tissue_barplot 
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_per_tissue_barplot.pdf", plot = l1_clust_prop_per_tissue_barplot,  width = 8, height = 8)

# G) Barplot of L1 cluster proportion per assigned sex at birth ----
all_samples_md <- l1_teaseq_obj@meta.data
asab_cell_tot <- all_samples_md %>%
  group_by(asab) %>%
  summarise(total = n(), .groups = "drop")
asab_cell_prop_df <- all_samples_md %>%
  group_by(l1_wnn_annot, asab) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(asab_cell_tot, by = "asab") %>%
  mutate(prop = count / total)
l1_clust_prop_asab_barplot <- ggplot(asab_cell_prop_df, aes(x = asab, y = prop, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = l1_clust_cols) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_asab_barplot 
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_asab_barplot.pdf", plot = l1_clust_prop_asab_barplot, width = 8, height = 8)

# H) Barplot of L1 cluster proportion per GEM ----
all_samples_md <- l1_teaseq_obj@meta.data
gem_cell_tot <- all_samples_md %>%
  group_by(GEM) %>%
  summarise(total = n(), .groups = "drop")
gem_cell_prop_df <- all_samples_md %>%
  group_by(l1_wnn_annot, GEM) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(gem_cell_tot, by = "GEM") %>%
  mutate(prop = count / total)
l1_clust_prop_per_gem_barplot <- ggplot(gem_cell_prop_df, aes(x = GEM, y = prop, fill = l1_wnn_annot)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = l1_clust_cols) +
  labs(y = "L1 Cluster Proportion",
       fill = "L1 Cluster") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 25, margin = margin(r = 10)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, 'sans', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, 'sans'),
    legend.text = element_text(size = 25, 'sans'),
    legend.title = element_text(size = 25, 'sans'),
    plot.margin = margin(t = 30, r = 5, b = 5, l = 5, unit = "pt")
  )
l1_clust_prop_per_gem_barplot 
ggsave("/filepath/fig1/fig1_supp/l1_wnn_clust_prop_barplots/l1_clust_prop_per_gem_barplot.pdf", plot = l1_clust_prop_per_gem_barplot, width = 10, height = 7)