# Sam Barnett Dubensky et al.
# Derek A. Oldridge & Laura A. Vella Labs at the Children's Hospital of Philadelphia
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned Tfh in humans
# Code and data visualization for Fig. S9 (related to Fig. 4)
# Fig. S9 - Derivation of spatially resolved GC Tfh signature and differential gene expression in FOXP3+ GC Th subset, related to Figure 4

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

# A) Combined spatial transcriptomic and TEAseq analyses resolve a core set of GC Tfh positional and lineage identity features ----

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

# B) FOXP3+ Tfh within the GC exhibit Treg-like transcriptional features and reduced GNG4 expression ----

# Import Xenium Prime L3 Tfh-only object with CN metadata
xp_tfh_obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l3_tfh_subset_obj_step10.rds')

# Inspect FOXP3 RNA expression distribution in L3 Tfh by violin plot
DefaultAssay(xp_tfh_obj) <- 'RNA'
VlnPlot(xp_tfh_obj, assay = 'RNA', features = 'FOXP3', sort = 'increasing', group.by = 'cn_gc_closed_annot') + 
  geom_hline(yintercept = 1) + NoLegend() + theme(axis.title = element_blank(), plot.title = element_blank())
VlnPlot(xp_tfh_obj, assay = 'RNA', features = 'FOXP3', group.by = 'cn_gc_closed_annot', pt.size = 0.01, sort = 'increasing') +
  geom_hline(yintercept = 1) + NoLegend() + theme(axis.title = element_blank(), plot.title = element_blank())

# Set threshold for positive vs negative FOXP3 expression at normalized expression value of 1
xp_tfh_obj$foxp3_rna_class <- ifelse(FetchData(xp_tfh_obj, vars = "rna_FOXP3") >= 1, "FOXP3_RNA_high", "FOXP3_RNA_low")
table(xp_tfh_obj$foxp3_rna_class) # Out of 235811 total cells in L3 Tfh object, 3288 (1.39%) are FOXP3+ vs 232523 (98.6%) FOXP3-

# Find DEGs between FOXP3+ vs FOXP3- Tfh that are spatially positioned within the GC (closed GC CN analysis, refer to Xenium Data Preprocessing Step 9)
xp_tfh_obj <- JoinLayers(xp_tfh_obj, assay = 'RNA')
xp_tfh_obj <- JoinLayers(xp_tfh_obj, assay = 'sketch')
Idents(xp_tfh_obj) <- 'foxp3_rna_class'
gc_closed_tfh_foxp3_hi_vs_lo_deg <- FindMarkers(subset(xp_tfh_obj, cn_gc_closed %in% '2'), # comparison subset on Tfh positioned within GCs
                                                ident.1 = 'FOXP3_RNA_high', 
                                                ident.2 = 'FOXP3_RNA_low', 
                                                assay = 'RNA', 
                                                min.pct = 0, logfc.threshold = 0) # removing filters for volcano plot visualization
gc_closed_tfh_foxp3_hi_vs_lo_deg$gene <- rownames(gc_closed_tfh_foxp3_hi_vs_lo_deg)

# Set threshold values and DEG color code for volcano plot - P < 1e-5 (0.00001) and log2FC > |0.25|
gc_closed_tfh_foxp3_hi_vs_lo_cols <- ifelse(
  gc_closed_tfh_foxp3_hi_vs_lo_deg$avg_log2FC < -0.25 & gc_closed_tfh_foxp3_hi_vs_lo_deg$p_val < 1e-5,  "#4b94d6",
  ifelse(
    gc_closed_tfh_foxp3_hi_vs_lo_deg$avg_log2FC >  0.25 & gc_closed_tfh_foxp3_hi_vs_lo_deg$p_val < 1e-5,  "#FE9A30",
    "lightgrey"                            
  )
)
names(gc_closed_tfh_foxp3_hi_vs_lo_cols)[gc_closed_tfh_foxp3_hi_vs_lo_cols == "lightgrey"] <- "NS"
names(gc_closed_tfh_foxp3_hi_vs_lo_cols)[gc_closed_tfh_foxp3_hi_vs_lo_cols == "#4b94d6"] <- "Down"
names(gc_closed_tfh_foxp3_hi_vs_lo_cols)[gc_closed_tfh_foxp3_hi_vs_lo_cols == "#FE9A30"] <- "Up in FOXP3+ GC Tfh"

# Find top FOXP3+ GC Tfh DEGs of interest to label
gc_closed_tfh_foxp3_hi_vs_lo_deg_pos <- gc_closed_tfh_foxp3_hi_vs_lo_deg %>% filter(avg_log2FC > 0.25) %>% filter(p_val < 1e-5) %>% arrange(desc(avg_log2FC)) %>% rownames()
gc_closed_tfh_foxp3_hi_vs_lo_deg_pos_labs <- gc_closed_tfh_foxp3_hi_vs_lo_deg_pos
gc_closed_tfh_foxp3_hi_vs_lo_deg_pos_labs <- gc_closed_tfh_foxp3_hi_vs_lo_deg_pos[ gc_closed_tfh_foxp3_hi_vs_lo_deg_pos != 'GBP5'] # GBP5 is above specified thresholds, but not included in labeled genes to visualize other genes of interest

# Find top FOXP3- GC Tfh DEGs of interest to label
gc_closed_tfh_foxp3_hi_vs_lo_deg_neg <- gc_closed_tfh_foxp3_hi_vs_lo_deg %>% filter(avg_log2FC < -0.25) %>% filter(p_val < 1e-5) %>% arrange(desc(avg_log2FC)) %>% rownames()
gc_closed_tfh_foxp3_hi_vs_lo_deg_neg_labs <- gc_closed_tfh_foxp3_hi_vs_lo_deg_neg

# Combine DEG lists for volcano plot
gc_closed_tfh_foxp3_hi_vs_lo_deg_vol_labs <- c(gc_closed_tfh_foxp3_hi_vs_lo_deg_pos_labs, gc_closed_tfh_foxp3_hi_vs_lo_deg_neg_labs)

# Assemble FOXP3+ vs FOXP3- GC Tfh DEG volcano plot
pdf('/filepath/fig4/fig4_supp/treg_spatial/gc_closed_tfh_foxp3_pos_vs_neg_deg_volcano.pdf', height = 6.5, width = 7)
EnhancedVolcano(gc_closed_tfh_foxp3_hi_vs_lo_deg,
                lab = rownames(gc_closed_tfh_foxp3_hi_vs_lo_deg),
                selectLab = gc_closed_tfh_foxp3_hi_vs_lo_deg_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'FOXP3+ vs FOXP3- GC Tfh DEGs',
                drawConnectors = FALSE,
                pCutoff = 1e-5,
                FCcutoff = 0.25,
                pointSize = 2.5,
                labSize = 4.5,
                colCustom= gc_closed_tfh_foxp3_hi_vs_lo_cols,
                legendPosition = 'none',
                labFace = 'italic',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  xlim(-1.7,4.69) # for visualization purposes, omitting features not meeting p-value threshold, as well as FOXP3 itself given the comparison
dev.off()

# Save complete DEG list for Data File S5
gc_closed_tfh_foxp3_hi_vs_lo_deg_export <- gc_closed_tfh_foxp3_hi_vs_lo_deg %>%
  rename(
    gene_symbol = gene,
    avg_log2fc_foxp3_pos_vs_neg_gctfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_foxp3_pos = pct.1,
    pct_pos_in_foxp3_neg = pct.2
  )
gc_closed_tfh_foxp3_hi_vs_lo_deg_export <- gc_closed_tfh_foxp3_hi_vs_lo_deg_export %>% arrange(p_val_adj)
gc_closed_tfh_foxp3_hi_vs_lo_deg_export <- gc_closed_tfh_foxp3_hi_vs_lo_deg_export %>% relocate(c('gene_symbol','avg_log2fc_foxp3_pos_vs_neg_gctfh','p_val_raw','p_val_adj'))
write_xlsx(gc_closed_tfh_foxp3_hi_vs_lo_deg_export, '/filepath/fig4/fig4_supp/treg_spatial/gc_closed_tfh_foxp3_hi_vs_lo_deg_export.xlsx')

# C-F) DAPI, CN, Treg/fr vs Tfh annotation, and FOXP3 spatial gene expression images ----

# Images exported from Xenium Explorer (v3.2.0)
