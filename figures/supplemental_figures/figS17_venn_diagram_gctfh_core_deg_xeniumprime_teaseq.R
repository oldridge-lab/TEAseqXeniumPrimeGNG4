# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S17 (related to Fig. 4M)
# Fig. S17 - Combined spatial transcriptomic and TEAseq analyses resolve a core set of GC Tfh positional and lineage identity features. 

# Set up R working environment ----

# Set working directory
setwd('/filepath/fig4/fig4_supp/')

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Matrix) # 1.7-0
library(tidyselect) # 1.2.1
library(Seurat) # 5.1.0
library(SeuratObject) # 5.0.2
library(data.table) # 1.15.4
library(readr) # 2.1.5
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4
library(stringr) # 1.5.1
library(presto) # 1.0.0
library(paletteer) # 1.6.0
library(RColorBrewer) # 1.1-3
library(EnhancedVolcano) # 1.22.0
library(tidyverse) # 2.0.0
library(BPCells) # 0.3.0
library(tibble) # 3.2.1
library(purrr) # 1.0.2
library(circlize) # 0.4.16
library(writexl) # 1.5.40
library(ggnewscale) # 0.5.0
library(readxl) # 1.4.3
library(openxlsx) # 4.2.8
library(scales) # 1.3.0
library(ggVennDiagram) # 1.5.2
library(pheatmap) # 1.0.12

# Import and filter DEG lists from three comparisons shown in Fig. 4M to define GC Tfh lineage and positioning features ----

# For each comparison, filtering on log2FC > 0 to select genes enriched in GC Tfh and Bonferroni-adjusted Wilcoxon rank-sum P-value < 0.05
# All comparisons use merged GC CN assignments (refer to Step 9 and Supplementary Materials)

# Xenium Prime - for all cells within closed GC CNs, what are DEGs enriched in Tfh vs all other cell types?
xp_gc_cn_gc_closed_tfh_vs_else_deg <- readRDS('/filepath/fig4/fig4m/xp_gc_cn_gc_closed_tfh_vs_else_deg.rds')
xp_gc_cn_gc_closed_tfh_vs_else_deg_filt <- xp_gc_cn_gc_closed_tfh_vs_else_deg %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) # XP - L2 Tfh versus all other L1/L2 cells within GC CN

# Xenium Prime - among L2 Tfh, what are DEGs between cells within closed GC vs nonGC CNs?
xp_tfh_gc_vs_nongc_cn_gc_closed_deg <- readRDS('/filepath/fig4/fig4m/xp_tfh_gc_vs_nongc_rna_deg_gc_closed_cn_l1_k40_n10.rds')
xp_tfh_gc_vs_nongc_cn_gc_closed_deg_filt <- xp_tfh_gc_vs_nongc_cn_gc_closed_deg %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) # XP - L2 Tfh in GC vs not in GC cellular neighborhoods

# TEAseq - what are DEGs enriched in L1 Tfh-like vs all other 14 L1 cell types? (RNA modality, restricted to 5101 genes measured in Xenium Prime experiment)
teaseq_tfh_vs_else_deg <- readRDS('/filepath/fig4/fig4m/teaseq_tfh_vs_else_deg.rds')
teaseq_tfh_vs_else_deg_filt <- teaseq_tfh_vs_else_deg %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) # TEAseq - L1 Tfh markers versus all other 14 subsets
teaseq_tfh_vs_else_deg_filt <- teaseq_tfh_vs_else_deg_filt %>% filter(gene %in% rownames(xp_tfh_gc_vs_nongc_cn_gc_closed_deg)) # Include only DEG within the 5101 gene XP panel
dim(teaseq_tfh_vs_else_deg_filt)

# Create Venn diagram comparing DEG lists from each comparison ----

# Assemble list of DEG sets for Venn diagram
gctfh_deg_venn_list <- list(Set1 = xp_tfh_gc_vs_nongc_cn_gc_closed_deg_filt$gene, # XP - Tfh in GCs versus Tfh positioned outside GCs
                            Set2 = xp_gc_cn_gc_closed_tfh_vs_else_deg_filt$gene, # XP - Tfh in GCs versus all other cell types in GCs
                            Set3 = teaseq_tfh_vs_else_deg_filt$gene) # L1 TEAseq object, RNA modality - Tfh-like cluster versus all other cell types

# Assemble Venn diagram
gctfh_deg_venn <- ggVennDiagram(gctfh_deg_venn_list, 
                                label_alpha = 0, 
                                label_size = 9,
                                category.names = c('GC Tfh vs nonGC Tfh','GC Tfh vs Other Cells in GC','Tfh vs Other Cell Types (TEAseq RNA)')) + # labels added in Illustrator
  scale_fill_gradient(low = "white", high = "#e49a78") +
  theme(legend.position = "none")
pdf('/filepath/fig4/fig4m/gctfh_deg_venn.pdf', height = 7, width = 7)
gctfh_deg_venn
dev.off()

# Find intersecting genes from all three comparisons and save list for Data File S5 ----
gctfh_common_genes <- Reduce(intersect, gctfh_deg_venn_list) 
gctfh_threeway_xp_teaseq_enriched_genes_df <- as.data.frame(gctfh_common_genes)
colnames(gctfh_threeway_xp_teaseq_enriched_genes_df) <- 'gctfh_threeway_xp_teaseq_enriched_genes'
write_xlsx(gctfh_threeway_xp_teaseq_enriched_genes_df, '/filepath/fig4/fig4_supp/xenium_gctfh_specificity_threeway_axis/gctfh_threeway_xp_teaseq_enriched_genes.xlsx')