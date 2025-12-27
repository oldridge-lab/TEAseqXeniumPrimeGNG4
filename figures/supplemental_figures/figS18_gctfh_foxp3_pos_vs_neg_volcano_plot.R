# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S18 (related to Fig. 4)
# Fig. S18 – FOXP3+ Tfh within the GC exhibit Treg-like transcriptional features and reduced GNG4 expression

# Set up R working environment

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