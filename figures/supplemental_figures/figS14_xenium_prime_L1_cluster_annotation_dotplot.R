# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S14 (related to Fig. 4H)
# Fig. S14 - Resolution of diverse immune and non-immune cell types in tonsil single-cell spatial transcriptomics.

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
library(pheatmap) # 1.0.12

# Import L1 Xenium Prime Seurat object
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l1_obj_step10.rds')

# RNA features for L1 annotation
xp_l1_annot_dotplot_features <- c('MS4A1','BANK1','CD72',
                                  'FCRL5','FCRL4','ITGAX',
                                  'MZB1','XBP1','SLAMF7',
                                  'CD38','LMO2','CD83','SERPINA9',
                                  'BCL6','TCF3','HMMR',
                                  'TYMS','RRM2','AURKB',
                                  'CD3E','CD4','TCF7','IL7R',
                                  'TOX2','GNG4','IL21','ICOS','PDCD1',
                                  'CD8A','CCL5','GZMA',
                                  'CLEC4C','LILRA4','IL3RA',
                                  'CPVL','AOAH','CLEC10A',
                                  'AIRE','CCR7','CD274','FSCN1',
                                  'CD68','MMP9','CTSL',
                                  'IL1B','PTGS2','TREM1',
                                  'HDC','MS4A2','KIT',
                                  'NOTCH3','PDGFRB','ANO1',
                                  'PECAM1','SELP','PLVAP',
                                  'GJA5','HEY1','SOX17',
                                  'PROX1','FOXC2','TIE1',
                                  'PI16','PDGFRA','COL5A1',
                                  'CCL19','CCL21','CXCL12',
                                  'S1PR3','CXCL13','GABRG3',
                                  'MSLN','CLDN18','VSTM2L',
                                  'CDH3','TP63','COL17A1',
                                  'FOXN1','DSG3','FAT2',
                                  'SPINK7','KLK6','CRNN')

# Specify cluster ordering for dotplot
xp_l1_annot_dotplot_order <- c('Naive B','MBC','ASC','LZ GCB','DZ GCB','Cycling',
                               'Naive T','nnCD4','Cytotox',
                               'PDC','DC','eTAC','Mono/Mac',
                               'Gran','Mast',
                               'Mural','Ven Endo','Art Endo','LEC',
                               'Trab Fib','FRC','FDC',
                               'Crypt Epi','Basal Epi','Int Epi','Surf Epi'
)
Idents(xp.obj) <- 'l1_annot'
xp.obj$l1_annot <- factor(
  xp.obj$l1_annot,
  levels = xp_l1_annot_dotplot_order
)

# Assemble dotplot
DefaultAssay(xp.obj) <- "RNA"
pdf("/filepath/fig4/fig4_supp/xp_l1_annot_dotplot_axes_top_right_up.pdf", width = 14, height = 6)
DotPlot(
  xp.obj,
  features = xp_l1_annot_dotplot_features,
  assay = "RNA",
  group.by = "l1_annot",
  cluster.idents = FALSE,
  cols = "PuOr",
  scale.by = "size",
  dot.scale = 4.5
) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, face = "italic",
      margin = margin(t = 2, b = 0), size = 9
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "right",
    legend.box.just = "center",
    legend.margin = margin(t = 10, b = -10),
    legend.spacing.x = unit(2, "cm"),
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
    legend.text  = element_text(size = 9, family = "sans", margin = margin(t = 3, b = 0)),
    legend.title = element_text(size = 9, family = "sans", margin = margin(t = 0, b = 3))
  ) +
  geom_vline(
    xintercept = seq(1.5, length(xp_l1_annot_dotplot_features) - 0.5, by = 1),
    color = "grey80", linetype = "solid", linewidth = 0.3
  ) +
  geom_hline(
    yintercept = seq(1.5, length(levels(xp.obj$l1_annot)) - 0.5, by = 1),
    color = "grey80",
    linewidth = 0.3
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      title = "Percent Expressing",
      title.hjust = 0.5
    ),
    color = guide_colourbar(
      order = 2,
      title.position = "top",
      title = "Scaled Expression",
      title.hjust = 0.5,
      barwidth = unit(6, "cm")
    )
  ) +
  scale_color_gradient2(
    low = RColorBrewer::brewer.pal(11, "PuOr")[10],
    mid = "white",
    high = RColorBrewer::brewer.pal(11, "PuOr")[2],
    midpoint = 0
  )
dev.off()