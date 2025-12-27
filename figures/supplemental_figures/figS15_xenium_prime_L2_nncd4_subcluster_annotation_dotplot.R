# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S15 (related to Fig. 4I)
# Fig. S15 - Level 2 subclustering analysis of nnCD4 T cells in tonsil Xenium Prime dataset reveals distinct Tfh versus nonTfh states. 

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

# Import Xenium Prime L2 nnCD4 T cell subclustering Seurat object
xp_nncd4t_obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l2_nncd4t_obj_step10.rds')

# Specify cluster ordering for dotplot
xp_l2_annot_dotplot_order <- c('Trm','Tconv','Tconv Myl','Treg/fr',
                               'Tcm 1','Tcm 2',
                               'Tfh Circ','Tfh CCR6','Tfh Mant',
                               'Tfh NFATC1','Tfh Myl','Tfh CXCL13','Tfh PRDM1',
                               'Tfh S1PR2','Tfh DZ','Tfh TOX2'
)
Idents(xp_nncd4t_obj) <- 'l2_nncd4t_annot'
xp_nncd4t_obj$l2_nncd4t_annot <- factor(
  xp_nncd4t_obj$l2_nncd4t_annot,
  levels = xp_l2_annot_dotplot_order
)

# Select RNA features to visualize for L2 nnCD4 T annotation dotplot
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_l2_annot_dotplot_features <- c('ITGAE','CXCR6','BHLHE40','RUNX3',
                                  'PRDM1','ID2',
                                  'IFNG','TBX21','CXCR3','ZEB2',
                                  'GATA3','CCR4','IL5','IL13',
                                  'RORA','RORC','IL17A','IL17B','IL22','IL26','CCR6','KLRB1','CCL20','PALLD','TNFRSF8','RUNX1','RUNX2',
                                  'FOXP3','IL2RA','CTLA4','TNFRSF1B','IL32','CYTOR','IL1R1','IKZF2','ENTPD1','PBXIP1','RTKN2',
                                  'CSF1R','CD68','ITGAX',
                                  'CCR7','SELL','GPR183','IL7R','TXNIP','KLF2','S1PR1',
                                  'CXCR5',
                                  'TYMS','PCNA','MCM4',
                                  'SKI',
                                  'UBASH3A','LAT','SEMA4D',
                                  'NFATC1','NR4A1','CD200','BTLA','EGR1','EGR2','IL4','IL21','CD40LG','PDCD1',
                                  'CXCL13',
                                  'CD38',
                                  'S1PR2','BCL6',
                                  'CXCR4',
                                  'TOX2','ASCL2','POU2AF1','POU3F1','B3GAT1','GNG4','XXYLT1','HEY1','MYB','CHGB','FZD3','CNIH3')

# Assemble L2 nnCD4 Annotation DotPlot
pdf("/filepath/fig4/fig4_supp/xp_l2_annot_dotplot_more_markers_legend_top.pdf", width = 14, height = 6)
DotPlot(
  xp_nncd4t_obj,
  features = xp_l2_annot_dotplot_features,
  assay = "RNA",
  group.by = "l2_nncd4t_annot",
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
    xintercept = seq(1.5, length(xp_l2_annot_dotplot_features) - 0.5, by = 1),
    color = "grey80", linetype = "solid", linewidth = 0.3
  ) +
  geom_hline(
    yintercept = seq(1.5, length(levels(xp_nncd4t_obj$l2_nncd4t_annot)) - 0.5, by = 1),
    color = "grey80", linewidth = 0.3
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      title = 'Percent Expressing',
      title.hjust = 0.5
    ),
    color = guide_colourbar(
      order = 2,
      title.position = "top",
      title = 'Scaled Expression',
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