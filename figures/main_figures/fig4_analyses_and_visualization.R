# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. 4
# Figure 4 - GNG4+ Tfh are primarily positioned within the germinal center light zone

# Set up R working environment ----

# Set working directory to Fig 4
setwd('/filepath/fig4/')

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
library(Rfast2) # 0.1.5.4
library(ggVennDiagram) # 1.5.2
library(ggrepel) # 0.9.5
library(openxlsx) # 4.2.8
library(scales) # 1.3.0

# Import Seurat objects for Xenium Prime and TEAseq datasets ----

# Xenium Prime object metadata includes annotations for L1-L3 cell type, cell neighborhood, and merged GC cell neighborhood
# All panels in Fig 4 relate to Xenium Prime spatial transcriptomics data - TEAseq L1 object is only used in Fig 4M

# Xenium Prime L1 object
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l1_obj_step10.rds')

# Xenium Prime L2 nnCD4 T subclustering object 
xp_nncd4t_obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l2_nncd4t_obj_step10.rds')

# Xenium Prime L3 Tfh-only subset object
xp_tfh_obj <- readRDS(file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l3_tfh_subset_obj_step10.rds')

# TEAseq L1 Seurat object (trimodal dimensionality reduction and L1 3WNN clustering of all tonsil and peripheral blood mononuclear cells, including Harmony integration across donors)
l1_teaseq_obj <- readRDS('/filepath/step13_bulk_harmony/l1_teaseq_obj.rds')

# Specify plot colors for visualization ----

# Specify colors for visualization of tonsil Xenium Prime Level 1 Seurat object (26 total clusters)
l1_xp_annot_colors <- c(
  "Naive B" = "#7b7f9c",
  "Naive T" = "#d699ac",
  "Mono/Mac" = "#c49980",
  "FRC" = "#90a6bf",
  "nnCD4" = "#c9744d", 
  "ASC" = "#b37cbf",
  "Int Epi" = "#e49a78",
  "LZ GCB" = "#f7f1e0", 
  "DZ GCB" = "#657ab0",
  "Cytotox" = "#dfcc78", 
  "Surf Epi" = "#d4b5e4", 
  "Cycling" = "#f4a6c7",
  "Basal Epi" = "#e1d4c7", 
  "Gran" = "#728782",
  "Ven Endo" = "#8ab1cc", 
  "Mural" = "#b6b7ba", 
  "Crypt Epi" = "#94a890",
  "DC" = "#9f5a6a",
  "PDC" = "#bad4f4",
  "MBC" = "#971a86", 
  "Art Endo" = "#d66d97",
  "FDC" = "#add5df",
  "LEC" = "#50664f",
  "Trab Fib" = "#e0c1b4", 
  "eTAC" = "#d19a11",
  "Mast" = '#69c9d1'
)

# Specify colors for visualization of tonsil Xenium Prime Level 2 Seurat object (16 total clusters)
nncd4t_cols_numb <- c(
  "0"  = "#117ab0",
  "1"  = "#f4a6c7",
  "2"  = "darkgrey",
  "3"  = "#c49980",
  "4"  = "#99c9b6",
  "5"  = "#c9744d",
  "6"  = "#bad4f4",
  "7"  = "#e49a78",
  "8"  = "#8ab1cc",
  "9"  = "#728782",
  "10" = "#e1d4c7",
  "11" = "#b37cbf",
  "12" = "#dfcc78",
  "13" = "#f7f1e0",
  "14" = "#94a890",
  "15" = "#e17794"
)
nncd4t_clust_cols <- setNames(nncd4t_cols_numb, l2_nncd4t_annot[names(nncd4t_cols_numb)])

# Set color palette for dot plots and scatter plots
puor_colors <- brewer.pal(11, "PuOr")

# Fig 4A - Xenium Prime experimental design schematic ----

# Prepared using BioRender and exported as PDF

# Fig 4B - Tonsil H&E Image ----

# H&E whole slide image, aligned to DAPI signal from Xenium Prime multimodal cell segmentation staining - performed in collaboration with Yutong Zhu from Oldridge Lab. Refer to code in 'step9_wsi_registration_cn_analysis' folder for details.
# TC653B used as representative donor (6YO M with SDB)
# Exported FOV from Xenium Explorer v3.2.0 and cropped in Adobe Illustrator v29.7.1

# Fig 4C - Tonsil DAPI Image ----

# TC653B used as representative donor (6YO M with SDB)
# Exported FOV from Xenium Explorer v3.2.0 and cropped in Adobe Illustrator v29.7.1

# Fig 4D - Tonsil BCL6 RNA Image ----

# TC653B used as representative donor (6YO M with SDB)
# Exported FOV from Xenium Explorer v3.2.0 and cropped in Adobe Illustrator v29.7.1
# Transcript Scaling = 12

# Fig 4E - Tonsil LMO2/AICDA RNA Image ----

# TC653B used as representative donor (6YO M with SDB)
# Exported FOV from Xenium Explorer v3.2.0 and cropped in Adobe Illustrator v29.7.1
# Transcript Scaling = 12

# Fig 4F - Tonsil GNG4 RNA Image ----

# TC653B used as representative donor (6YO M with SDB)
# Exported FOV from Xenium Explorer v3.2.0 and cropped in Adobe Illustrator v29.7.1
# Transcript Scaling = 12

# Fig 4G - Tonsil GNG4/AICDA RNA Image ----

# TC653B used as representative donor (6YO M with SDB)
# Magnified FOV of B-F, exported from Xenium Explorer v3.2.0 and cropped in Adobe Illustrator v29.7.1
# Transcript Scaling = 12

# Fig 4H - L1 RNA UMAP of Tonsil XP Clusters ----

# Level 1 annotations for 26 clusters from FindMarkers resolution 1.0
xp_l1_annot <- c(
  '0'  = 'Naive B', # Relative to MBC, less ITGAX 
  '1'  = 'Naive T', # Mixed naive CD4 and CD8 T cell features
  '2'  = 'Mono/Mac', # Mixed monocyte and mature macrophage-related genes
  '3'  = 'FRC', # Fibroblastic reticular cells
  '4'  = 'nnCD4', # non-naive CD4 T cells, Tfh-like signature including TOX2, BCL6, and GNG4
  '5'  = 'ASC', # Antibody-secreting cell phenotype
  '6'  = 'Int Epi', # Found to be spatially positioned between surface and basal epithelium including within crypts using Xenium Explorer
  '7'  = 'LZ GCB', # relative to DZ-like cluster, more CXCR5, LMO2, FCER2 (CD23), CD83, CD40, CD86, CD72, SEMA7A, BCL2A1, EGR3, CIITA
  '8'  = 'DZ GCB', # relative to LZ-like cluster, more CXCR4, AICDA, HMMR, TCF3, POLH, LIG4, and proliferation-related genes
  '9'  = 'Cytotox', # Including features of NK cells and CD8 T cell
  '10' = 'Surf Epi', # Found to be spatially positioned at periphery of epithelium including within crypts using Xenium Explorer
  '11' = 'Cycling', # Largely proliferation-related genes, but also features of several lineages, especially B and T cells
  '12' = 'Basal Epi', # Found to be spatially positioned at basal layer of epithelium including within crypts using Xenium Explorer
  '13' = 'Gran', # Granulocytes, particularly neutrophil-related genes in contrast to c25 Mast-like cells
  '14' = 'Ven Endo', # Venous-like endothelium
  '15' = 'Mural', # Mural-like cells found to be positioned around or near endothelial cells using Xenium Explorer
  '16' = 'Crypt Epi', # Epithelial cell cluster primarily found in crypts but not exterior epithelium using Xenium Explorer 
  '17' = 'DC', # Dendritic cells. Relative to c24 eTAC-like cells, more CPVL ITGAX AOAH SAMHD1 ITGB2 CLEC10A CLEC4A CYBB FOS AHNAK NAIP ANPEP SCIMP CD4 FCGR2A CLEC7A ANXA1 MYCL CIITA
  '18' = 'PDC', # Plasmacytoid dendritic cells
  '19' = 'MBC', # Memory B cells
  '20' = 'Art Endo', # Arterial-like endothelial cells
  '21' = 'FDC', # Follicular dendritic cells
  '22' = 'LEC', # Lymphatic endothelial-like cells
  '23' = 'Trab Fib', # Fibroblasts enriched in trabecular positioning i.e. thick connective tissue bands running through tonsil samples
  '24' = 'eTAC', # Extrathymic AIRE-expressing antigen-presenting cells. Relative to c17 DC cells, more AIRE AMP3 FSCN1 CD83 RGS1 CDKN1A IL7R ADAM12 LAD1 CCL22 DUSP5 KIF2A NRXN2 LY75.
  '25' = 'Mast' # Separate granulocyte cluster from 'Gran' with mast cell signature
)

# For data and rationale supporting cluster annotations, refer to Fig. S14 and Data File S4 for key DEGs, Fig. S16 for spatial distribution, and Supplementary Materials for references and discussion

# Apply L1 annotations to full object
DefaultAssay(xp.obj) <- 'RNA'
xp.obj$l1_annot <- xp.obj$full_bulk_clust_1.0_res
Idents(xp.obj) <- 'l1_annot'
xp.obj <- RenameIdents(xp.obj, xp_l1_annot)
xp.obj$l1_annot <- Idents(xp.obj)

# Annotated L1 XP UMAP - rasterized due to very high cell number
pdf('/filepath/fig4/fig4h/xp_l1_umap_annot_repel.pdf')
Idents(xp.obj) <- 'l1_annot'
DefaultAssay(xp.obj) <- 'RNA'
DimPlot(xp.obj, group.by = 'l1_annot', cols = l1_xp_annot_colors, alpha = 0.2, label = TRUE, label.box = TRUE, label.size = 4, repel = TRUE, raster = TRUE) + 
  NoLegend() + coord_fixed() + NoAxes() + theme(plot.title = element_blank())
dev.off()

# Fig 4I - L2 nnCD4 T Cell Subclustering UMAP ----

# Annotate L2 nnCD4 T cell subclusters
l2_nncd4t_annot <- c(
  '0'  = 'Tfh DZ',
  '1'  = 'Tfh Mant',
  '2'  = 'Tcm 1',
  '3'  = 'Tfh CXCL13',
  '4'  = 'Tfh TOX2',
  '5'  = 'Tfh S1PR2',
  '6'  = 'Treg/fr',
  '7'  = 'Tconv Myl', 
  '8'  = 'Tfh CCR6',
  '9'  = 'Trm', 
  '10' = 'Tconv', 
  '11' = 'Tcm 2',
  '12' = 'Tfh Myl',
  '13' = 'Tfh NFATC1',
  '14' = 'Tfh Circ', 
  '15' = 'Tfh PRDM1'
)

# Refer to Fig. 4J, Fig. S15, and Data File S4 for key DEGs, Fig. S16 for spatial distribution, and Supplementary Materials for references and additional discussion of rationale behind cluster annotations.

# UMAP for XP L2 nnCD4 T cell subclusters with offset cluster labels
Idents(xp_nncd4t_obj) <- 'l2_nncd4t_annot'
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
nncd4t_umap_dimplot <- DimPlot(xp_nncd4t_obj, reduction = 'full.umap.sketch.1', alpha = 0.7, raster = TRUE, 
                               group.by = 'l2_nncd4t_annot', cols = nncd4t_clust_cols) + 
                               coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
nncd4t_umap_coords <- FetchData(xp_nncd4t_obj, vars = c("fullumapsketch1_1", "fullumapsketch1_2")) %>%
  mutate(l2_nncd4t_annot = Idents(xp_nncd4t_obj),
         fill_col = nncd4t_clust_cols[as.character(l2_nncd4t_annot)])
nncd4t_umap_df <- nncd4t_umap_coords %>% 
  group_by(l2_nncd4t_annot, fill_col) %>% 
  summarize(UMAP_1 = median(fullumapsketch1_1), UMAP_2 = median(fullumapsketch1_2), .groups = 'drop')
nncd4t_umap_lab <- nncd4t_umap_dimplot +
  new_scale_color() + 
  geom_label_repel(
    data = nncd4t_umap_df,
    aes(x = UMAP_1, y = UMAP_2, label = l2_nncd4t_annot, fill = fill_col, color = 'black'),
    size = 4.25,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines"),
    force = 1,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  NoLegend()
pdf('/filepath/fig4/fig4i/l2_nncd4_umap_raster_repel_labs.pdf', width = 6, height = 4)
nncd4t_umap_lab
dev.off()

# Set color palette for FeaturePlots below
or_cols <- brewer.pal(9, "Oranges")

# Using downsampled cells in 'sketch' assay for FeaturePlots to improve visualization, rather than rasterizing all cells from 'RNA' assay
DefaultAssay(xp_nncd4t_obj) <- 'sketch'

# CXCR5
pdf('/filepath/fig4/fig4i/nncd4t_umap_cxcr5.pdf', width = 6, height = 5)
FeaturePlot(xp_nncd4t_obj, raster = FALSE, order = TRUE, slot = 'data', reduction = 'umap.sketch',
            features = "CXCR5", cols = c(or_cols[1], or_cols[8])) +
  guides(
    color = guide_colorbar(
      barwidth  = unit(0.5, "cm"),  
      barheight = unit(6.5, "cm")
    )
  ) +
  theme(
    legend.position      = "left",
    legend.box.margin    = margin(t = 0, r = -3, b = 0, l = 0, unit = "cm"),
    legend.box.spacing   = unit(0, "cm"),
    plot.title = element_blank()
  ) +
  NoAxes() + 
  coord_fixed()
dev.off()

# PDCD1
pdf('/filepath/fig4/fig4i/nncd4t_umap_pdcd1.pdf', width = 6, height = 5)
FeaturePlot(xp_nncd4t_obj, raster = FALSE, order = TRUE, slot = 'data', reduction = 'umap.sketch',
            features = "PDCD1", cols = c(or_cols[1], or_cols[8])) +
  guides(
    color = guide_colorbar(
      barwidth  = unit(0.5, "cm"),  
      barheight = unit(6.5, "cm")
    )
  ) +
  theme(
    legend.position      = "left",
    legend.box.margin    = margin(t = 0, r = -3, b = 0, l = 0, unit = "cm"),
    legend.box.spacing   = unit(0, "cm"),
    plot.title = element_blank()
  ) +
  NoAxes() + 
  coord_fixed()
dev.off()

# BCL6
pdf('/filepath/fig4/fig4i/nncd4t_umap_bcl6.pdf', width = 6, height = 5)
FeaturePlot(xp_nncd4t_obj, raster = FALSE, order = TRUE, slot = 'data', reduction = 'umap.sketch',
            features = "BCL6", cols = c(or_cols[1], or_cols[8])) +
  guides(
    color = guide_colorbar(
      barwidth  = unit(0.5, "cm"),  
      barheight = unit(6.5, "cm")
    )
  ) +
  theme(
    legend.position      = "left",
    legend.box.margin    = margin(t = 0, r = -3, b = 0, l = 0, unit = "cm"),
    legend.box.spacing   = unit(0, "cm"),
    plot.title = element_blank()
  ) +
  NoAxes() + 
  coord_fixed()
dev.off()

# GNG4
pdf('/filepath/fig4/fig4i/nncd4t_umap_gng4.pdf', width = 6, height = 5)
FeaturePlot(xp_nncd4t_obj, raster = FALSE, order = TRUE, slot = 'data', reduction = 'umap.sketch',
            features = "GNG4", cols = c(or_cols[1], or_cols[8])) +
  guides(
    color = guide_colorbar(
      barwidth  = unit(0.5, "cm"),  
      barheight = unit(6.5, "cm")
    )
  ) +
  theme(
    legend.position      = "left",
    legend.box.margin    = margin(t = 0, r = -3, b = 0, l = 0, unit = "cm"),
    legend.box.spacing   = unit(0, "cm"),
    plot.title = element_blank()
  ) +
  NoAxes() + 
  coord_fixed()
dev.off()

# CCR6
pdf('/filepath/fig4/fig4i/nncd4t_umap_ccr6.pdf', width = 6, height = 5)
FeaturePlot(xp_nncd4t_obj, raster = FALSE, order = TRUE, slot = 'data', reduction = 'umap.sketch',
            features = "CCR6", cols = c(or_cols[1], or_cols[8])) +
  guides(
    color = guide_colorbar(
      barwidth  = unit(0.5, "cm"),  
      barheight = unit(6.5, "cm")
    )
  ) +
  theme(
    legend.position      = "left",
    legend.box.margin    = margin(t = 0, r = -3, b = 0, l = 0, unit = "cm"),
    legend.box.spacing   = unit(0, "cm"),
    plot.title = element_blank()
  ) +
  NoAxes() + 
  coord_fixed()
dev.off()

# Fig 4J - Annotation DotPlot of Tonsil XP L2 nnCD4 T Cell Clusters ----

# Select DEGs for dotplot visualization
xp_l2_annot_dotplot_features <- c('ITGAE','CXCR6',
                                  'PRDM1','ID2','IFNG',
                                  'FOXP3','IL2RA','CXCR3',
                                  'CCR7','KLF2','GPR183','IL7R',
                                  'CXCR5','CCR6',
                                  'LAT','FYN','BTLA','CD40LG',
                                  'NFATC1','NR4A1','IL21','IL4','PDCD1',
                                  'CXCL13',
                                  'CD38',
                                  'S1PR2','BCL6',
                                  'CXCR4',
                                  'TOX2','B3GAT1','GNG4')


# Set cluster order for dotplot visualization
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

# Assemble L2 nnCD4 Annotation DotPlot
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_l2_annot_dotplot_axes_top <- DotPlot(xp_nncd4t_obj, features = xp_l2_annot_dotplot_features, assay = 'RNA', 
                                        group.by = 'l2_nncd4t_annot', cluster.idents = FALSE, 
                                        cols = 'PuOr', scale.by = 'size', dot.scale = 4.5) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', margin = margin(t = 2, b = 0), size = 9),
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
    legend.text = element_text(size = 9,  family = 'sans', margin = margin(t = 2, b = 0)),
    legend.title = element_text(size = 9, family = 'sans', margin = margin(t = 0, b = 3))
  ) +
  geom_vline(xintercept = seq(1.5, length(xp_l2_annot_dotplot_features) - 0.5, by = 1),
             color     = "grey80",
             linetype  = "solid",
             linewidth      = 0.3) +
  geom_hline(
    yintercept = seq(1.5, length(levels(xp_nncd4t_obj$l2_nncd4t_annot)) - 0.5, by = 1),
    color      = "grey80",
    size       = 0.3
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      title = '', # Percent Expressed label added in Illustrator
      title.hjust = 0.5,
      label.position = "top",
      label.theme = element_text(margin = margin(t = 0, b = -0.5))
    ),
    color = guide_colourbar(
      order = 2,
      title.position = "top",
      label.position = "top",
      title = '', # Average Expression label added in Illustrator
      title.hjust = 0.5,
      barwidth = unit(2.5, "cm"),
      barheight = unit(0.3, "cm"),
      label.theme = element_text(margin = margin(t = -0.4, b = 3))
    )
  ) +
  scale_color_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  )
ggsave("/filepath/fig4/fig4j/xp_l2_annot_dotplot_axes_top.pdf",
       plot = xp_l2_annot_dotplot_axes_top, width = 6, height = 5.75)

# Dotplot further annotated in Illustrator to emphasize Tfh versus nonTfh groups

# Fig 4K - Xenium Explorer image of cell neighborhoods in tonsil ----

# Cell neighborhood (CN) analysis used L1 cluster annotations as input (not L2 T cell subclusters), k = 40 neighbors, and n = 10 neighborhoods
# Original method from Schürch, Christian M. et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. Cell, Volume 182, Issue 5, 1341 - 1359.e19
# Code adapted from related GitHub page - https://github.com/nolanlab/NeighborhoodCoordination
# CN generation for tonsil Xenium Prime dataset was performed in collaboration with Yutong Zhu from Oldridge Lab.
# For details regarding CN generation refer to Python notebooks in the Xenium Data Processing Step 9 folder.
# CNs were annotated as detailed in Xenium Data Processing Step 10 and Fig. S16 code.

# Panel 4K shows a representative image of cells in tonsil colored by assigned CN using sample TC653B.

# Fig 4L -  Positive predictive value of each gene for GC spatial positioning in Tfh, ranked from lowest to highest ----

# Get RNA expression matrix per cell in L3 Tfh-only subset object
tfh_rna_expr_mtx <- GetAssayData(xp_tfh_obj, assay = DefaultAssay(xp_tfh_obj), slot = "data")

# Get cell barcodes for L3 Tfh spatially positioned within the closed GC CN
gctfh_barcodes <- xp_tfh_obj$cn_gc_closed_annot == "CN2 GC"

# Define threshold for positive gene expression as 1
pos_thresh  <- tfh_rna_expr_mtx > 1

# For all 5101 genes in the XP panel, get total number of Tfh with positive expression
n_pos <- Matrix::rowSums(pos_thresh)

# For each gene, get number of Tfh with positive expression that are also spatially positioned within GC
n_pos_gc <- Matrix::rowSums(pos_thresh[, gctfh_barcodes, drop = FALSE])

# Find percentage of Tfh expressing each gene that are also spatially positioned within GC
pct_in_gc <- 100 * n_pos_gc / n_pos

# Assemble data table for GC Tfh gene positive predictive value analysis
gctfh_gene_ppv <- tibble(
  gene = rownames(tfh_rna_expr_mtx),
  n_pos = as.integer(n_pos),
  n_pos_in_GC = as.integer(n_pos_gc),
  pct_pos_in_GC = as.numeric(pct_in_gc)
)

# Rank all 5101 genes measured from highest to lowest PPV for Tfh spatial positioning in GC
gctfh_gene_ppv <- gctfh_gene_ppv %>% arrange(pct_pos_in_GC) %>% mutate(gene_ord = factor(gene, levels = gene))

# Rename table columns and save for Data File S5
xp_tfh_rna_gc_pct_all_genes <- gctfh_gene_ppv %>%
  rename(
    gene_rna = gene,
    n_tfh_gene_pos = n_pos,
    n_tfh_gene_pos_in_gc = n_pos_in_GC,
    pct_tfh_gene_pos_in_gc = pct_pos_in_GC
  ) %>%
  arrange(desc(pct_tfh_gene_pos_in_gc)) # Rank all 5101 genes measured from highest to lowest PPV for Tfh spatial positioning in GC
write.xlsx(xp_tfh_rna_gc_pct_all_genes, '/filepath/fig4/fig4l/xp_tfh_rna_gc_pct_all_genes.xlsx')

# Select genes to visualize in ranked gene plot
label_genes <- c(
  "GNG4","TOX2","B3GAT1","IL21","BCL6","PDCD1","CXCR5","CCR6","GPR183",
  "KLF2","RORC","FOXP3","IL4","IL10","S1PR2","IL17A",
  "CXCR3","CCR4","TBX21","GATA3","IFNG","IL5","IL13","CXCL13","ICOS",
  "TIGIT","IL7R","THY1","MAF","IRF4"
)

# Assemble dataframe for plot labels to emphasize GNG4
lab_df <- gctfh_gene_ppv %>%
  filter(gene %in% label_genes) %>%
  mutate(
    is_gng4 = gene == "GNG4",
    label_size = ifelse(is_gng4, 8, 5), # make GNG4 label size for visualization
    label_face = ifelse(is_gng4, "bold.italic", "italic")
  )

# Create ranked gene PPV plot
gctfh_gene_ppv_plot <- ggplot(gctfh_gene_ppv, aes(gene_ord, pct_pos_in_GC)) +
  geom_point(size = 0.15) +
  geom_hline(yintercept = 50, color = "black", linewidth = 0.75, linetype = 'dashed') +
  ggrepel::geom_text_repel(
    data = lab_df,
    aes(
      label = gene,
      color = is_gng4,
      size = label_size,
      fontface = label_face
    ),
    max.overlaps  = Inf,
    box.padding   = 0.4,
    point.padding = 0.6,
    nudge_y = 1,
    force = 3,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0 
  ) +
  scale_color_manual(
    values = c(`TRUE` = "#FF6464", `FALSE` = "black"),
    guide = "none"
  ) +
  scale_size_identity(guide = "none") +
  scale_y_continuous(
    labels = scales::label_number(),  
    breaks = seq(0, 100, by = 10),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "% Tfh Expressing Marker Gene that are within GC") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = 'sans', size = 15),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save plot
ggsave(
  "/filepath/fig4/fig4l/pct_tfh_in_gc_by_gene_ranking.pdf",
  gctfh_gene_ppv_plot, width = 5.5, height = 5, units = "in"
)

# Next, we directly compared PPV for GNG4 versus panel of other genes conventionally associated with GC Tfh
# Here, PPV is calculated for each gene on a per-donor basis - in contrast to calculations for data shown in 4L, where PPV was assessed using all Tfh across samples regardless of sample origin 
# Statistics are shown in Data File S5

# Get list of donors, in stable order for functions below
donor_list <- sort(unique(xp_tfh_obj$donor))

# Get vector of donor identifiers for all L3 Tfh cells
donor_vec <- xp_tfh_obj$donor

# For each donor, determine percentage of Tfh expressing a given gene that localize to the GC CN
pct_by_donor <- do.call(cbind, lapply(donor_list, function(dd) {
  
  # Get cells for given donor
  ii  <- donor_vec == dd
  
  # Define positive expression for each gene as greater than 1 in each cell
  pos <- tfh_rna_expr_mtx[, ii, drop = FALSE] > 1
  
  # Get percentages for each gene in each donor
  100 * Matrix::rowSums(pos[, gctfh_barcodes[ii], drop = FALSE]) / Matrix::rowSums(pos)
  
}))

# Assign gene names to output matrix
rownames(pct_by_donor) <- rownames(tfh_rna_expr_mtx)

# Assign donor names to output matrix
colnames(pct_by_donor) <- donor_list

# Genes for comparison with GNG4
gng4_comp_panel <- c(
  "BCL6","TOX2","B3GAT1","S1PR2","IL4","ASCL2","POU2AF1",
  "CXCL13","CD200","PDCD1","TIGIT","IL21",
  "ICOS","CXCR5","SH2D1A","CD40LG","CXCR4","SLAMF6","BTLA",
  "CHGB","TOX","MAF","NFATC1"
)

# Filter RNA expression matrix from above to genes of interest
gng4_comp_panel <- gng4_comp_panel[gng4_comp_panel %in% rownames(tfh_rna_expr_mtx)]

# Extract per donor percentages for GNG4 (gene 'x' in function below)
x <- pct_by_donor["GNG4", ]

# Function to compare percent GC-localized Tfh for each gene in panel versus GNG4 on a per-donor basis
stats_df <- do.call(rbind, lapply(gng4_comp_panel, function(g) {
  
  y <- pct_by_donor[g, ] # Extract per-donor percentages for other gene 'g'
  
  # Perform paired t-test between GNG4 and other gene
  tt <- t.test(x, y,
               paired = TRUE, # within donor
               alternative = "greater")  # GNG4 > gene
  
  # Perform paired Wilcoxon signed-rank test between GNG4 and other gene
  pw <- wilcox.test(x, y,
                    paired = TRUE,
                    alternative = "greater",
                    exact = FALSE)
  
  # Save summary statistics and p-values for each gene compared
  data.frame(
    gene = g,
    mean_pct_GNG4 = mean(x, na.rm = TRUE),
    mean_pct_gene = mean(y, na.rm = TRUE),
    mean_diff = mean(x - y, na.rm = TRUE), # Difference in percentage between GNG4 and given gene
    p_t_paired = tt$p.value,
    p_wilcox_pair = pw$p.value,
    stringsAsFactors = FALSE
  )
}))

# Correct for multiple comparisons using Benjamini-Hochberg FDR method
stats_df$FDR_t <- p.adjust(stats_df$p_t_paired, method = "BH")
stats_df$FDR_wilcox <- p.adjust(stats_df$p_wilcox_pair, method = "BH")

# Rename columns and export for Data File S5
xp_tfh_rna_gc_pct_stats <- stats_df %>% rename(
  other_gene = gene,
  GNG4_mean_pct = mean_pct_GNG4,
  other_gene_mean_pct = mean_pct_gene,
  mean_diff_pct_gng4_vs_other = mean_diff,
  p_val_raw_t_paired = p_t_paired,
  FDR_val_t_paired = FDR_t,
  p_val_raw_wilcox_paired = p_wilcox_pair,
  FDR_val_wilcox_paired = FDR_wilcox) %>%
  arrange(desc(other_gene_mean_pct)) %>%
  relocate(c('GNG4_mean_pct','other_gene','other_gene_mean_pct','mean_diff_pct_gng4_vs_other','p_val_raw_t_paired','FDR_val_t_paired','p_val_raw_wilcox_paired','FDR_val_wilcox_paired'))
write.xlsx(xp_tfh_rna_gc_pct_stats, '/filepath/fig4/fig4l/xp_tfh_rna_gc_pct_stats.xlsx')

# Fig 4M - Three-way scatterplot of DEGs enriched in GC Tfh based on positioning and cell type comparisons ----

# Join assay layers for DEG analysis
xp.obj <- JoinLayers(xp.obj, assay = 'RNA')
xp.obj <- JoinLayers(xp.obj, assay = 'sketch')
l1_teaseq_obj <- JoinLayers(l1_teaseq_obj, assay = 'RNA')

# Create groups for Tfh versus nonTfh subclusters
tfh_clust_names <- l2_nncd4t_annot[startsWith(l2_nncd4t_annot, "Tfh")]
nontfh_clust_names <- l2_nncd4t_annot[!startsWith(l2_nncd4t_annot, "Tfh")]

# XP - for all cells within closed GC CNs, what are DEGs enriched in Tfh vs all other cell types?
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_nncd4t_annot'
xp_gc_cn_gc_closed_tfh_vs_else_deg <- FindMarkers(subset(xp.obj, subset = cn_gc_closed %in% c("2")), # numerical assignment for closed GC CN
                                                  assay = 'RNA',
                                                  ident.1 = tfh_clust_names, # all Tfh clusters resolved in L2 subclustering analysis
                                                  ident.2 = c('Art Endo','ASC','Basal Epi','Crypt Epi','Cycling','Cytotox','DC','DZ GCB','eTAC','FDC','FRC','Gran','Int Epi','LEC','LZ GCB','Mast','MBC','Mono/Mac','Mural','Naive B','Naive T','PDC','Surf Epi','Trab Fib','Ven Endo',
                                                              nontfh_clust_names),
                                                  only.pos = FALSE,
                                                  logfc.threshold = 0,
                                                  min.pct = 0)
xp_gc_cn_gc_closed_tfh_vs_else_deg$gene <- rownames(xp_gc_cn_gc_closed_tfh_vs_else_deg) # 5101 genes measured in Xenium Prime experiment
saveRDS(xp_gc_cn_gc_closed_tfh_vs_else_deg, '/filepath/fig4/fig4m/xp_gc_cn_gc_closed_tfh_vs_else_deg.rds')

# XP - among L2 Tfh, what are DEGs between cells within closed GC vs nonGC CNs?
DefaultAssay(xp_tfh_obj) <- 'RNA'
Idents(xp_tfh_obj) <- 'cn_gc_closed'
xp_tfh_gc_vs_nongc_cn_gc_closed_deg <- FindMarkers(
  xp_tfh_obj,
  ident.1 = c('2'), # numerical assignment for closed GC CN
  ident.2 = c('0','1','3','4','6','7','8','9'), # all other CNs
  only.pos = FALSE,
  min.pct = 0, logfc.threshold = 0,
  assay = 'RNA')
xp_tfh_gc_vs_nongc_cn_gc_closed_deg$gene <- rownames(xp_tfh_gc_vs_nongc_cn_gc_closed_deg) # 5101 genes measured in Xenium Prime experiment
saveRDS(xp_tfh_gc_vs_nongc_cn_gc_closed_deg, '/filepath/fig4/fig4m/xp_tfh_gc_vs_nongc_rna_deg_gc_closed_cn_l1_k40_n10.rds')

# TEAseq - what are DEGs enriched in L1 Tfh-like vs all other 14 L1 cell types? (RNA modality, restricted to 5101 genes measured in Xenium Prime experiment)
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
teaseq_tfh_vs_else_deg <- FindMarkers(l1_teaseq_obj,
                                      assay = 'RNA',
                                      ident.1 = 'Tfh-like',
                                      ident.2 = ,
                                      only.pos = FALSE,
                                      logfc.threshold = 0,
                                      min.pct = 0)
teaseq_tfh_vs_else_deg$gene <- rownames(teaseq_tfh_vs_else_deg) # 36601 genes in whole-transcriptome TEAseq RNA assay
saveRDS(teaseq_tfh_vs_else_deg, '/filepath/fig4/fig4m/teaseq_tfh_vs_else_deg.rds')

# Renaming columns for each output DEG list above

# XP - for all cells within closed GC CNs, what are DEGs enriched in Tfh vs all other cell types?
xp_gc_cn_gc_closed_tfh_vs_else_deg <- readRDS('/filepath/fig4/fig4m/xp_gc_cn_gc_closed_tfh_vs_else_deg.rds')
colnames(xp_gc_cn_gc_closed_tfh_vs_else_deg)
colnames(xp_gc_cn_gc_closed_tfh_vs_else_deg) <- c('pval_tfh_spec','fc_tfh_spec','pct_tfh','pct_nontfh','pvaladj_tfhspec','gene')
dim(xp_gc_cn_gc_closed_tfh_vs_else_deg) # 5101 genes measured in Xenium Prime experiment

# XP - among L2 Tfh, what are DEGs between cells within closed GC vs nonGC CNs?
xp_tfh_gc_vs_nongc_cn_gc_closed_deg <- readRDS('/filepath/fig4/fig4m/xp_tfh_gc_vs_nongc_rna_deg_gc_closed_cn_l1_k40_n10.rds')
colnames(xp_tfh_gc_vs_nongc_cn_gc_closed_deg)
colnames(xp_tfh_gc_vs_nongc_cn_gc_closed_deg) <- c('pval_gctfh','fc_gctfh','pct_gctfh','pct_nongctfh','pvaladj_gctfh','gene')
dim(xp_tfh_gc_vs_nongc_cn_gc_closed_deg) # 5101 genes measured in Xenium Prime experiment

# In the L1 TEAseq dataset, what are DEG enriched in Tfh-like cells vs all other cell types in tonsil and PBMC?
teaseq_tfh_vs_else_deg <- readRDS('/filepath/fig4/fig4m/teaseq_tfh_vs_else_deg.rds')
colnames(teaseq_tfh_vs_else_deg)
colnames(teaseq_tfh_vs_else_deg) <- c('pval_tfh_teaseq','fc_tfh_teaseq','pct_tfh_teaseq','pct_nontfh_teaseq','pvaladj_tfh_teaseq','gene')
dim(teaseq_tfh_vs_else_deg)  # 36601 genes in whole-transcriptome TEAseq RNA assay
teaseq_tfh_vs_else_deg_xp_feats_only <- teaseq_tfh_vs_else_deg %>% filter(gene %in% rownames(xp_gc_cn_gc_closed_tfh_vs_else_deg)) # limit comparison to genes measured by Xenium Prime assay
dim(teaseq_tfh_vs_else_deg_xp_feats_only) # now restricted to 5101 genes measured in Xenium Prime experiment

# Assemble three-way enrichment dataframe
gctfh_enr_df <- merge(xp_tfh_gc_vs_nongc_cn_gc_closed_deg, xp_gc_cn_gc_closed_tfh_vs_else_deg, by = 'gene')
gctfh_enr_df <- merge(gctfh_enr_df, teaseq_tfh_vs_else_deg_xp_feats_only, by = 'gene')
dim(gctfh_enr_df) # 5101 genes

# Filtering plotted DEGs based on Bonferroni-adjusted Wilcoxon rank-sum test P < 0.05
gctfh_enr_df_filt <- gctfh_enr_df %>% 
  filter(
    pvaladj_gctfh      < 0.05,
    pvaladj_tfhspec    < 0.05,
    pvaladj_tfh_teaseq < 0.05
  )

# DEGs to label on enrichment plot
label_genes <- c(
  "CCR6", "S1PR1", "PLAC8", "GPR183", "P2RX1", "PDGFD", "P2RY12", "INHBB",
  "PEG10", "GRIK4", "WNK2", "DAB1", "HAL", "B3GAT1", "CHGB", "TOX2",
  "KSR2", "PDCD1", "LIF", "IL21", "CTLA4", "KLRB1", "LAG3", "IL7R", "F5",
  "GIMAP5", "CCR7", "LINC00861",
  'TIGIT','IL4'
)

# Build label data frame for enrichment plot, in which all genes are italicized, and GNG4 is bold
lab_df <- gctfh_enr_df_filt %>%
  filter(gene %in% c(label_genes, "GNG4")) %>%
  mutate(
    label_expr = ifelse(
      gene == "GNG4", "bold(italic('GNG4'))", # GNG4 bold and italic
      paste0("italic('", gene, "')") # All other genes italic only
    ),
    label_size = ifelse(gene == "GNG4", 6.5, 4) # Larger GNG4 text size
  )

# Assemble three-way GC Tfh enrichment plot
gctfh_enr_scatter_plot_filt <- ggplot(
  gctfh_enr_df_filt,
  aes(x = fc_gctfh, y = fc_tfh_spec)
  ) +
  
  # Lines at origin
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  
  # Point layer for X,Y coordinates of GC vs nonGC Tfh DEGs and GC Tfh vs Other GC Cell DEGs
  geom_point(
    aes(fill = fc_tfh_teaseq), # Points colored by Tfh vs nonTfh in L1 TEAseq dataset RNA modality
    shape = 21,
    color  = "black",
    stroke = 0.2,
    size = 2.2
  ) +
  
  # Label layer for genes
  geom_label_repel(
  data = lab_df,
  aes(
    label = label_expr,
    size  = label_size,
    fill  = fc_tfh_teaseq
  ),
  parse = TRUE, # Parsing italic and bold gene names set above
  color = "black",
  label.size = 0.25,
  box.padding = 0.75,
  point.padding = 0.5,
  force = 1,
  max.overlaps = Inf,
  segment.color = "grey40",
  min.segment.length = 0.1,
  show.legend = FALSE) +
  scale_size_identity() + 
  scale_fill_gradient2(
    low = puor_colors[10],
    mid = "white",
    high = puor_colors[2],
    midpoint = 0
  ) +
  
  # Legend for L1 TEAseq fold-change
  guides(
    fill = guide_colorbar(
      title = NULL, # added manually in Illustrator
      label = FALSE,
      barwidth = 8,
      barheight = 0.8
    )
  ) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right"
  )

# Save plot
pdf('/filepath/fig4/fig4m/gctfh_enr_scatter_plot_filt.pdf', width = 7, height = 5.25)
gctfh_enr_scatter_plot_filt
dev.off()

# Export DEG lists for each comparison visualized on plot as single spreadsheet for Data File S5

# XP - Tfh in GC versus Tfh not in GC
xp_tfh_gc_vs_nongc_cn_gc_closed_deg <- readRDS('/filepath/fig4/fig4m/xp_tfh_gc_vs_nongc_rna_deg_gc_closed_cn_l1_k40_n10.rds')
xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers <- xp_tfh_gc_vs_nongc_cn_gc_closed_deg %>%
  rename(
    avg_log2fc_gctfh_vs_nongctfh = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_gctfh = pct.1,
    pct_pos_in_nongctfh = pct.2
  )
xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers <- xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers %>% arrange(desc(avg_log2fc_gctfh_vs_nongctfh), p_val_adj)
xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers <- xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers %>% relocate(c('gene_symbol','avg_log2fc_gctfh_vs_nongctfh','p_val_raw','p_val_adj'))
saveRDS(xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers, '/filepath/fig4/fig4m/xp_tfh_gc_vs_nongc_cn_gc_closed_deg_markers.rds')

# XP - Tfh in GC vs Other Cells in GC
xp_gc_cn_gc_closed_tfh_vs_else_deg <- readRDS('/filepath/fig4/fig4m/xp_gc_cn_gc_closed_tfh_vs_else_deg.rds')
xp_gc_cn_gc_closed_tfh_vs_else_deg_markers <- xp_gc_cn_gc_closed_tfh_vs_else_deg %>%
  rename(
    avg_log2fc_tfh_vs_others_in_gc = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_gctfh = pct.1,
    pct_pos_in_other_gc_cells = pct.2
  )
xp_gc_cn_gc_closed_tfh_vs_else_deg_markers <- xp_gc_cn_gc_closed_tfh_vs_else_deg_markers %>% arrange(desc(avg_log2fc_tfh_vs_others_in_gc), p_val_adj)
xp_gc_cn_gc_closed_tfh_vs_else_deg_markers <- xp_gc_cn_gc_closed_tfh_vs_else_deg_markers %>% relocate(c('gene_symbol','avg_log2fc_tfh_vs_others_in_gc','p_val_raw','p_val_adj'))
saveRDS(xp_gc_cn_gc_closed_tfh_vs_else_deg_markers, '/filepath/fig4/fig4m/xp_gc_cn_gc_closed_tfh_vs_else_deg_markers.rds')

# TEAseq RNA - L1 Tfh-like cells vs all other L1 cell types
teaseq_tfh_vs_else_deg <- readRDS('/filepath/fig4/fig4m/teaseq_tfh_vs_else_deg.rds')
teaseq_tfh_vs_else_deg_markers <- teaseq_tfh_vs_else_deg %>%
  rename(
    avg_log2fc_tfhlike_vs_others = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_tfh = pct.1,
    pct_pos_in_others = pct.2
  )
teaseq_tfh_vs_else_deg_markers <- teaseq_tfh_vs_else_deg_markers %>% arrange(desc(avg_log2fc_tfhlike_vs_others), p_val_adj)
teaseq_tfh_vs_else_deg_markers <- teaseq_tfh_vs_else_deg_markers %>% relocate(c('gene_symbol','avg_log2fc_tfhlike_vs_others','p_val_raw','p_val_adj'))
saveRDS(teaseq_tfh_vs_else_deg_markers, '/filepath/fig4/fig4m/teaseq_tfh_vs_else_deg_markers.rds')

# Export DEG lists as single spreadsheet for Data File S5
clust_marker_df_dir   <- '/filepath/fig4/fig4m/'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/fig4/fig4m/xenium_gctfh_specificity_threeway_axis_deg.xlsx")  

# Fig 4N - Spatial positioning of all Tfh versus GNG4+ Tfh versus LZ and DZ GC B cells ----

# Panel shows spatial positioning of all L3 Tfh, GNG4+ L3 Tfh cells only, and L1 LZ or DZ GC B cells
# Cell annotations were visualized using Xenium Explorer in a magnified FOV of Panel 4K, as indicated by dashed box
# Image was further annotated in Adobe Illustrator

# In code below we define GNG4+ vs GNG4- Tfh cells in the Xenium Prime dataset, and extract cell annotations to import into Xenium Explorer

# Inspect GNG4 RNA expression distribution in L3 Tfh by violin plot
VlnPlot(xp_tfh_obj, assay = 'RNA', features = 'GNG4', sort = 'increasing') + 
  geom_hline(yintercept = 1) + NoLegend() + theme(axis.title = element_blank(), plot.title = element_blank())

# Set threshold for positive vs negative GNG4 expression at normalized expression value of 1
xp_tfh_obj$gng4_rna_class <- ifelse(FetchData(xp_tfh_obj, vars = "rna_GNG4") >= 1, "GNG4_RNA_high", "GNG4_RNA_low")
table(xp_tfh_obj$gng4_rna_class) # 52532 (22.3%) pos vs 183279 (77.7%) neg (out of 235811 total)

# Get 'GNG4_RNA_high' versus 'GNG4_RNA_low' character vector with L3 Tfh cell barcode names
gctfh_gng4_annot_vec <- as.character(xp_tfh_obj@meta.data$gng4_rna_class) # GNG4 class
names(gctfh_gng4_annot_vec) <- rownames(xp_tfh_obj@meta.data) # L3 Tfh cell barcodes

# Get L1 object metadata and combined L1/L2 object annotations
l1_obj_md  <- xp.obj@meta.data
l1_l2_labs <- as.character(l1_obj_md$l1_l2_nncd4t_annot)

# Create new annotation vector in which original L2 Tfh labels are replaced with GNG4 class, while remaining L2 nnCD4 T cell and L1 cell labels are retained
l1_l2_annot_with_tfh_gng4_class <- l1_l2_labs # combined L1 and L2 nnCD4 T cell annotations
names(l1_l2_annot_with_tfh_gng4_class) <- rownames(l1_obj_md) # all cell barcodes
l1_l2_annot_with_tfh_gng4_class[names(gctfh_gng4_annot_vec)] <- gctfh_gng4_annot_vec # GNG4 class for Tfh cells only

# Add combined L1/L2 and L3 Tfh GNG4 class annotations to L1 object
xp.obj$l1_l2_annot_with_tfh_gng4_class <- l1_l2_annot_with_tfh_gng4_class
table(xp.obj@meta.data$l1_l2_annot_with_tfh_gng4_class)

# TC653B - extract annotations and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'l1_l2_annot_with_tfh_gng4_class'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/fig4/fig4n/tc653b_l1_with_l2_tfh_gng4_class_annotations.csv",
  row.names = FALSE,
  quote     = FALSE
)

# Fig 4O - Smoothed line plot of gene expression in Tfh as a function of distance from the nearest LZ GC B cell ----

# In code below, we extract spatial coordinates and gene expression for each Tfh cell, determine the distance of each cell from the nearest L1 LZ GC B cell, and create a smoothed line plot of average gene expression versus distance

# Extract RNA expression from Tfh subset object for genes of interest
DefaultAssay(xp_tfh_obj) <- 'RNA' # 235811 cells
tfh_rna_expr_df <- FetchData(xp_tfh_obj, vars = c('GNG4','BCL6','TOX2','B3GAT1','S1PR2','THY1','POU2AF1','ASCL2','TIGIT','XXYLT1','CD200','CXCR5','PDCD1','IL21','ICOS','CCR6','GPR183','RORA','CCR7','KLF2','RORC','KLRB1','IL7R','IL17A'))
saveRDS(tfh_rna_expr_df, '/filepath/fig4/fig4o/tfh_rna_expr_df.rds')

# Get centroid coordinates for all segmented cell types in all samples
coords <- map_dfr(
  Images(xp.obj),
  ~ GetTissueCoordinates(xp.obj, image = .x, scale = "hires") %>%
    mutate(image = .x)
)

# Rename Tfh coordinate matrix row and column names
tfh_coords <- coords %>%
  filter(cell %in% colnames(xp_tfh_obj)) %>%
  rename(imagecol = x, imagerow = y)

# Get LZ GC B cell coordinate matrix
lzgcb_subset <- subset(xp.obj, subset = full_bulk_clust_1.0_res == "7") # numerical cluster for LZ GC B cells in L1 object
lzgcb_coords <- coords %>%
  filter(cell %in% colnames(lzgcb_subset)) %>%
  rename(imagecol = x, imagerow = y)

# For each donor sample, get Euclidean distance between each Tfh and nearest LZ GC B cell
dist_df <- map_dfr(
  unique(tfh_coords$image), # will loop over all six samples
  function(fov) { # one fov per sample
    cd4_sub <- filter(tfh_coords, image == fov) # L2 Tfh cells in sample
    gcb_sub <- filter(lzgcb_coords, image == fov) # L1 LZ GCB cells in sample
    
    # Find single nearest LZ GC B cell for each Tfh
    knn_out <- get.knnx(
      data  = as.matrix(gcb_sub %>% select(imagecol, imagerow)),
      query = as.matrix(cd4_sub %>% select(imagecol, imagerow)),
      k     = 1
    )
    
    # Create data table with one row per Tfh cell and distance from the nearest LZ GCB cell
    tibble( 
      cell = cd4_sub$cell,
      dist_to_GCB = knn_out$nn.dist[,1]
    )
  }
)

# Join LZ GCB distance and Tfh RNA expression dataframes
tfh_rna_expr_df_dist <- tfh_rna_expr_df %>%
  rownames_to_column("cell") %>%
  left_join(tfh_coords, by = "cell") %>%
  left_join(dist_df, by = "cell")

# Select genes of interest and visualization colors for expression versus distance line plot
genes_of_interest <- c("GNG4","BCL6","CXCR5","PDCD1","IL21",'RORC','CCR6','IL7R','GPR183')
gene_colors <- c(
  "CXCR5"  = "black",
  "PDCD1"  = "black",
  "BCL6"   = "black",
  "IL21"   = "black",
  "RORC"   = "grey",
  "CCR6"   = "grey",
  "GPR183" = "grey",
  "IL7R"   = "grey",
  "GNG4"   = "#FF6464"
)

# Define kernel smoothing function to create smoothed line plot of average Tfh gene expression versus distance
# Original function from Derek A. Oldridge

' Kernel Smoothing Function
#'
#' Performs kernel smoothing on input data (xi, yi) to produce
#' smoothed estimates (yo) at specified output points (xo).
#'
#' @param xi Numeric vector of input x-values (independent variable).
#' @param yi Numeric vector of input y-values (dependent variable).
#' @param xo Numeric vector of output x-values where smoothing is applied.
#' @param kernel Character string specifying kernel type:
#'        - "move_average": assigns equal weight within a window of width const_tuning.
#'        - "gaussian": assigns weights following a Gaussian distribution.
#' @param const_tuning Numeric bandwidth parameter that can be tuned to control smoothness.
#'        Larger values give smoother results, smaller values retain more detail.
#'
#' @return A data frame with:
#'         - xo: output x-values
#'         - yo: smoothed y-values

kernel_smooth <- function(xi, yi, xo, kernel = c("move_average", "gaussian"), const_tuning) {
  
  # Choose the kernel type and define its weighting function
  if(kernel == "move_average") {
    # Moving average kernel: equal weight if within window of width const_tuning, otherwise 0
    weight_func <- function(xr, x, width = const_tuning) {
      abs(x - xr) <= width / 2
    }
  }
  if(kernel == "gaussian") {
    # Gaussian kernel: weights decrease exponentially with squared distance
    weight_func <- function(xr, x, sigma = const_tuning) {
      exp(-(x - xr)^2 / (2 * sigma^2))
    }
  }
  
  # Apply smoothing for each output point (xo)
  yo <- sapply(
    xo, # loop over each output point
    function(xr) {
      # Compute weights for all xi relative to current xo = xr
      w <- weight_func(xr, xi, const_tuning); 
      # Return weighted average of yi values
      return(sum(w*yi) / sum(w))
    }
  )
  
  # Return smoothed values as a data frame
  return(data.frame(xo, yo))
}

# Define parameters to run kernel smoothing function on gene expression versus distance data
xo <- seq(0, 100, by = 0.1)  # Smoothing grid (um)
bandwidth_um <- 5 # Kernel sigma
kernel_type <- "gaussian"

# Pivot expression versus distance dataframe from wide to long format and filter
tfh_rna_expr_df_dist_long <- tfh_rna_expr_df_dist %>%
  # Given that the expected radius of a lymphocyte is ~ 3 microns, remove cells less than 6 microns from nearest LC GC B cell due to likely segmentation error (https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=15&id=111438)
  dplyr::filter(dist_to_GCB > 6) %>% 
  dplyr::select(cell, dist_to_GCB, all_of(genes_of_interest)) %>%
  tidyr::pivot_longer(
    cols = all_of(genes_of_interest),
    names_to = "gene",
    values_to = "expression"
  )

# Execute smoothing function 
smooth_df <- tfh_rna_expr_df_dist_long %>%
  filter(!is.na(dist_to_GCB), !is.na(expression),
         dist_to_GCB >= 6, dist_to_GCB <= 100) %>%  # restricting smoothing grid to 6-100 microns
  group_by(gene) %>%
  reframe({
    out <- kernel_smooth(
      xi = dist_to_GCB,
      yi = expression,
      xo = seq(6, 100, by = 0.1), # restricting smoothing grid to 6-100 microns
      kernel = "gaussian",
      const_tuning = bandwidth_um
    )
    tibble(dist_um = out$xo, expr_smoothed = out$yo)
  })

# Define layering order for genes on line plot
smooth_df$gene <- factor(
  smooth_df$gene,
  levels = c("RORC","CCR6","IL7R","GPR183","BCL6","CXCR5","PDCD1","IL21","GNG4")
)
smooth_df <- dplyr::arrange(smooth_df, gene != "GNG4") # GNG4 on top layer

# Assemble smoothed line plot
pdf('/filepath/fig4/fig4o/tfh_expr_vs_lzgcb_dist_smooth.pdf', width = 5.3, height = 6.25)
ggplot(smooth_df, aes(x = dist_um, y = expr_smoothed, color = gene)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = gene_colors, name = NULL, limits = names(gene_colors)) +
  scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
  coord_cartesian(xlim = c(6, 100)) +
  scale_y_continuous(expand = c(0.01, 0)) + 
  labs(x = "Distance to nearest LZ GCB (µm)", y = "Smoothed RNA expression (raw data, 6–100 µm)") +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 13),
    legend.text      = element_text(size = 13),
    legend.position  = "right",
    legend.key       = element_blank(),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 10),
    axis.title       = element_blank(),
    axis.text        = element_text(size = 13),
    axis.ticks       = element_line(size = 0.5),
    axis.ticks.length= unit(0.2, "cm")
  ) +
  NoLegend()
dev.off()

# Fig 4P - DEGs between GNG4+ vs GNG4- Tfh spatially positioned within GCs ----

# Computing DEGs between GNG4+ vs GNG4- L3 Tfh, subset on cells assigned to merged GC CNs
xp_tfh_obj <- JoinLayers(xp_tfh_obj, assay = 'RNA')
DefaultAssay(xp_tfh_obj) <- 'RNA'
Idents(xp_tfh_obj) <- 'gng4_rna_class'
gctfh_gng4_hi_vs_lo_deg <- FindMarkers(subset(xp_tfh_obj, cn_gc_closed %in% c(2)), # comparison subset on Tfh positioned within GCs
                                       ident.1 = 'GNG4_RNA_high', 
                                       ident.2 = 'GNG4_RNA_low', 
                                       assay = 'RNA', 
                                       min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) # removing filters for volcano plot visualization
gctfh_gng4_hi_vs_lo_deg$gene <- rownames(gctfh_gng4_hi_vs_lo_deg)

# Set threshold values and DEG color code for volcano plot - P < 1e-10 (0.0000000001) and log2FC > |0.25|
gctfh_gng4_hi_vs_lo_cols <- ifelse(
  gctfh_gng4_hi_vs_lo_deg$avg_log2FC < -0.25 & gctfh_gng4_hi_vs_lo_deg$p_val < 1e-10,  "#4b94d6", 
  ifelse(
    gctfh_gng4_hi_vs_lo_deg$avg_log2FC >  0.25 & gctfh_gng4_hi_vs_lo_deg$p_val < 1e-10,  "#FF6464",
    "lightgrey"                            
  )
)
names(gctfh_gng4_hi_vs_lo_cols)[gctfh_gng4_hi_vs_lo_cols == "lightgrey"] <- "NS"
names(gctfh_gng4_hi_vs_lo_cols)[gctfh_gng4_hi_vs_lo_cols == "#4b94d6"] <- "Down"
names(gctfh_gng4_hi_vs_lo_cols)[gctfh_gng4_hi_vs_lo_cols == "#FF6464"] <- "Up in GNG4+ GC Tfh"

# Find top GNG4+ GC Tfh DEGs of interest to label
gctfh_gng4_hi_vs_lo_deg_pos <- gctfh_gng4_hi_vs_lo_deg %>% filter(avg_log2FC > 0.25) %>% filter(p_val < 1e-10) %>% arrange(desc(avg_log2FC)) %>% rownames()
gctfh_gng4_hi_vs_lo_deg_pos_labs <- gctfh_gng4_hi_vs_lo_deg_pos
gctfh_gng4_hi_vs_lo_deg_pos_labs <- gctfh_gng4_hi_vs_lo_deg_pos_labs[ 
  !gctfh_gng4_hi_vs_lo_deg_pos_labs %in% c('PTPN14', 'GRID1', 'PTHLH','XXYLT1','FZD3') # all of these DEGs are above specified thresholds, but not included in labeled genes to better visualize other genes of interest
]

# Find top GNG4- GC Tfh DEGs of interest to label
gctfh_gng4_hi_vs_lo_deg_neg <- gctfh_gng4_hi_vs_lo_deg %>% filter(avg_log2FC < -0.25) %>% filter(p_val < 1e-10) %>% arrange(desc(avg_log2FC)) %>% rownames()
gctfh_gng4_hi_vs_lo_deg_neg_labs <- gctfh_gng4_hi_vs_lo_deg_neg
gctfh_gng4_hi_vs_lo_deg_neg_labs <- gctfh_gng4_hi_vs_lo_deg_neg_labs[ 
  !gctfh_gng4_hi_vs_lo_deg_neg_labs %in% c('BMPER', 'AOAH') # all of these DEGs are above specified thresholds, but not included in labeled genes to better visualize other genes of interest
]

# Combine DEG lists for volcano plot
gctfh_gng4_hi_vs_lo_deg_vol_labs <- c(gctfh_gng4_hi_vs_lo_deg_pos_labs, gctfh_gng4_hi_vs_lo_deg_neg_labs)

# Assemble volcano plot
pdf('/filepath/fig4/fig4p/gctfh_gng4_pos_vs_neg_deg_volcano.pdf', height = 7.25, width = 7)
EnhancedVolcano(gctfh_gng4_hi_vs_lo_deg,
                lab = rownames(gctfh_gng4_hi_vs_lo_deg),
                selectLab = gctfh_gng4_hi_vs_lo_deg_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'GNG4+ vs GNG4- GC Tfh DEGs',
                drawConnectors = FALSE,
                pCutoff = 1e-10,
                FCcutoff = 0.25,
                pointSize = 2.5,
                labSize = 5,
                colCustom= gctfh_gng4_hi_vs_lo_cols,
                legendPosition = 'none',
                labFace = 'italic',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf,
                axisLabSize = 14
                ) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        axis.line = element_line(linewidth = 0.45),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  xlim(-1.88,1.56) # for visualization purposes, omitting features not meeting p-value threshold, as well as GNG4 itself given the comparison
dev.off()

# Save complete DEG list for Data File S5
gctfh_gng4_hi_vs_lo_deg_save <- gctfh_gng4_hi_vs_lo_deg %>%
  rename(
    gene_symbol = gene,
    avg_log2fc_gng4_pos_vs_neg_gctfh = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_gng4_pos = pct.1,
    pct_pos_in_gng4_neg = pct.2
  )
gctfh_gng4_hi_vs_lo_deg_save <- gctfh_gng4_hi_vs_lo_deg_save %>% arrange(p_val_adj)
gctfh_gng4_hi_vs_lo_deg_save <- gctfh_gng4_hi_vs_lo_deg_save %>% relocate(c('gene_symbol','avg_log2fc_gng4_pos_vs_neg_gctfh','p_val_raw','p_val_adj'))
write_xlsx(gctfh_gng4_hi_vs_lo_deg_save, '/filepath/fig4/fig4p/gctfh_gng4_hi_vs_lo_deg_save.xlsx')