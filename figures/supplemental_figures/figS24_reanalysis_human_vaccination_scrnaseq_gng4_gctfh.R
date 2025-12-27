# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. S24 (related to Fig. 5I & 5J)
# Fig. S24 - GNG4 expression is specifically induced in human GC-like Tfh during both SARS-CoV-2 and influenza vaccine responses.

# Set up R working environment ----

# Set working directory
setwd('/filepath/fig5/fig5_supp/')

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
library(scales) # 1.3.0
library(rlang) # 1.1.4

# Specify color palette for dot plots
rdbu_colors <- brewer.pal(11, "RdBu")

# Specify colors for numbered clusters 0-26 in reanalysis of Borcherding et al. SARS2 study
cols_borcherding <- c(
  "0"  = "#add5df",
  "1"  = "#8ab1cc",
  "2"  = "#94a890",
  "3"  = "#d66d97",
  "4"  = "#e49a78",
  "5"  = "#657ab0",
  "6"  = "#c49980",
  "7"  = "#72d1b4",
  "8"  = "#dfcc78",
  "9"  = "#b6b7ba",
  "10" = "#d4b5e4",
  "11" = "#69c9d1",
  "12" = "#e1d4c7",
  "13" = "#c9a3c1",
  "14" = "#f4a6c7",
  "15" = "#e17794",
  "16" = "#dfcc78",
  "17" = "#728782",
  "18" = "#c49980",
  "19" = "#b6b7ba",
  "20" = "#d699ac",
  "21" = "#50664f",
  "22" = "#971a86",
  "23" = "#69c9d1",
  "24" = "#f7f1e0",
  "25" = "#c9744d",
  "26" = "#b37cbf"
)

# Specify colors for numbered clusters 0-26 in reanalysis of Schattgen et al. IIV study
cols_schattgen <- c(
  "0"  = "#8ab1cc",
  "1"  = "#add5df",
  "2"  = "#94a890",
  "3"  = "#c9a3c1",
  "4"  = "#e49a78",
  "5"  = "#657ab0",
  "6"  = "#bad4f4",
  "7"  = "#72d1b4",
  "8"  = "#c9744d",
  "9"  = "#d66d97",
  "10" = "#50664f",
  "11" = "#d4b5e4",
  "12" = "#e1d4c7",
  "13" = "#69c9d1",
  "14" = "#e17794",
  "15" = "#f4a6c7",
  "16" = "#dfcc78",
  "17" = "#728782",
  "18" = "#c49980",
  "19" = "#b6b7ba",
  "20" = "#d699ac",
  "21" = "#50664f",
  "22" = "#971a86",
  "23" = "#69c9d1",
  "24" = "#f7f1e0",
  "25" = "#c9744d",
  "26" = "#b37cbf"
)

# Specify colors for Level 2 Tfh subclustering object in reanalysis of Schattgen et al. IIV study
tfh_type_cols <- c(
  'GC' = '#e17794',
  'IL10 TFH' = '#657ab0',
  'Treg' = '#e49a78',
  'cycling' = '#94a890',
  'pre/memory' = '#dfcc78',
  'naive' = '#add5df'
)

# Import Seurat objects for reanalysis ----

# For Fig S24A-S24C and Fig 5I -
# Reanalysis of Borcherding et al. CD4+ T cells exhibit distinct transcriptional phenotypes in the lymph nodes and blood following mRNA vaccination in humans. Nat Immunol. 2024. PMID 39164479. https://pubmed.ncbi.nlm.nih.gov/39164479/
# Import T cell subclustering Seurat object from SARS2 mRNA vaccination LN FNA scRNAseq study
# Seurat object .rds file obtained from Figure2 tab of CellPilot data repository (https://cellpilot.emed.wustl.edu/)
covid_tcell_obj <- readRDS('/filepath/reanalysis/borcherding_mrna_vax/BorcherdingFig2.rds')

# For Fig S24D-S24K and Fig 5J -
# Reanalysis of Schattgen et al. Influenza vaccination stimulates maturation of the human T follicular helper cell response. Nat Immunol. 2024. PMID 39164477. https://pubmed.ncbi.nlm.nih.gov/39164477/
# Import Seurat objects for T cell clustering and Level 2 Tfh subclustering from IIV LN FNA scRNAseq study 
# Seurat object .rds files ('intergrated_Tcells_harmony_bydonor.rds' and 'intergrated_Tfh_harmony_bydonor.rds') obtained from Zenodo data repository (Version 3, Oct 20 2023, https://zenodo.org/records/12611325)
flu_tcell_obj <- readRDS('/filepath/reanalysis/intergrated_Tcells_harmony_bydonor.rds')
flu_tfh_obj <- readRDS('/filepath/reanalysis/intergrated_Tfh_harmony_bydonor.rds')

# Fig S24A - UMAP of CD4 and CD8 T cell clustering in Borcherding et al. SARS2 study ----

Idents(covid_tcell_obj) <- 'seurat_clusters' # original clusters 0-11, with c3 annotated as CD4+ GC Tfh
pdf('/filepath/fig5/fig5_supp/covid_tcell_clust_numb_labs_dimplot.pdf', width = 6, height = 6)
DimPlot(covid_tcell_obj, label = TRUE, order = TRUE, label.box = TRUE, cols = cols_borcherding, pt.size = 1) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()

# c0-c11 cluster annotations derive from original study, listed to right of UMAP figure in Illustrator

# Fig S24B - UMAP of GNG4 RNA expression in T cells (reanalysis of Borcherding et al. study) ----

# UMAP without legend
Idents(covid_tcell_obj) <- 'seurat_clusters'
pdf('/filepath/fig5/fig5_supp/covid_tcell_clust_numb_labs_gng4_feat_plot.pdf', width = 6, height = 6)
FeaturePlot(covid_tcell_obj, features = 'GNG4', pt.size = 1, label = TRUE, order = TRUE, cols = c('grey','coral2')) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()

# UMAP with legend, cropped to legend and placed in Fig S24B UMAP figure panel using Illustrator
Idents(covid_tcell_obj) <- 'seurat_clusters'
pdf('/filepath/fig5/fig5_supp/covid_tcell_clust_numb_labs_gng4_feat_plot_legend.pdf', width = 6, height = 6)
FeaturePlot(covid_tcell_obj, features = 'GNG4', pt.size = 1, label = TRUE, order = TRUE, cols = c('grey','coral2')) + coord_fixed() + theme(plot.title = element_blank())
dev.off()

# Fig S24C - DotPlot of RNA expression for GNG4 and other DEGs across T cell clusters (reanalysis of Borcherding et al. study) ----

pdf('/filepath/fig5/fig5_supp/covid_tcell_seurat_clust_dotplot.pdf', height = 5, width = 9)
Idents(covid_tcell_obj) <- 'seurat_clusters'
DotPlot(covid_tcell_obj, features = c('CD4','CD8A','MKI67','LEF1','CCR7','IL7R','KLF2','GPR183','FOXP3','IL2RA','IL10','CXCR5','PDCD1','TIGIT','BCL6','TOX2','GNG4'), dot.scale = 5.25,
        group.by = 'seurat_clusters', cluster.idents = TRUE) +
  scale_color_gradient2(
    low = rdbu_colors[11],
    mid = "white",
    high = rdbu_colors[2],
    midpoint = 0    
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'italic', family = 'sans', angle = 45, hjust = 1))
dev.off()

# Fig S24D-S24K - Mapping L2 GC Tfh cells to L1 T cell clustering object and exploring specificity of GNG4 RNA expression (reanalysis of Schattgen et al. study) ----

# Note - as provided, the Level 1 T cell object does not have a specifically annotated 'GC Tfh' cluster, whereas the Level 2 Tfh subclustering object does ('GC')
# Therefore, in steps below we examine expression of GNG4 and other GC-like Tfh features in both the L2 (S24D-F) and L1 objects (S24G-I), then map the 'GC' cluster from L2 to L1 (S24J-K)

# Using Schattgen et al. annotations, label each cell in L2 Tfh subclustering object as GC Tfh 'yes' or 'no'
table(flu_tfh_obj$Tfh_type) # annotations from authors
flu_tfh_obj$gctfh_annot <- flu_tfh_obj$Tfh_type
Idents(flu_tfh_obj) <- 'gctfh_annot'
gctfh_annot <- c(
  'GC' = 'yes',
  'IL10 TFH' = 'no',
  'Treg' = 'no',
  'cycling' = 'no',
  'pre/memory' = 'no',
  'naive' = 'no'
)
flu_tfh_obj <- RenameIdents(flu_tfh_obj, gctfh_annot)
flu_tfh_obj$gctfh_annot <- Idents(flu_tfh_obj)
table(flu_tfh_obj$gctfh_annot)

# Retrieve cell barcodes for 'GC' Tfh vs all other annotated cells in L2 Tfh subset object
gctfh_annot_vec <- as.character(flu_tfh_obj@meta.data$gctfh_annot)
names(gctfh_annot_vec) <- rownames(flu_tfh_obj@meta.data)

# Apply L2 GCTfh vs nonGCTfh/nonTfh labels to L1 T cell object
flu_tcell_obj <- AddMetaData(flu_tcell_obj, metadata = gctfh_annot_vec, col.name = "gctfh_annot")
table(flu_tcell_obj$gctfh_annot, useNA = "ifany") # NA cells were not included in Tfh subset object and thus we group them into nonTfh annotation
flu_tcell_obj@meta.data$gctfh_annot[is.na(flu_tcell_obj@meta.data$gctfh_annot)] <- "no"
table(flu_tcell_obj$gctfh_annot, useNA = "ifany") # NA cells now relabeled 'no' in GC Tfh status - 2142 GCTfh now annotated in L1 object, as in L2 object

# Make new metadata column for L1 T cell object with 'GC Tfh' cells annotated, and original cluster numbers for all other cells retained ('seurat_clusters')
flu_tcell_obj$gctfh_vs_others_annot <- ifelse(
  flu_tcell_obj$gctfh_annot == "yes",
  "GC Tfh",
  as.character(flu_tcell_obj$seurat_clusters)
)

# Fig S24D - UMAP of Level 2 Tfh subclustering object from Schattgen et al. study ----

Idents(flu_tfh_obj) <- 'Tfh_type' # annotations from authors
table(flu_tfh_obj$Tfh_type)
pdf('/filepath/fig5/fig5_supp/flu_tfh_dimplot.pdf', width = 6, height = 6)
DimPlot(flu_tfh_obj, group.by = 'Tfh_type', cols = tfh_type_cols, order = TRUE, pt.size = 2, label.box = TRUE, label = TRUE) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()
  
# Fig S24E - UMAP of GNG4 RNA expression in L2 Tfh subclustering object (reanalysis of Schattgen et al. study) ----

# UMAP without legend
pdf('/filepath/fig5/fig5_supp/flu_tfh_gng4_feat_plot_nolegend.pdf', width = 6, height = 6)
FeaturePlot(flu_tfh_obj, features = 'GNG4', pt.size = 2, order = TRUE, cols = c('grey','coral2')) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()

# UMAP with legend, then cropped to legend and placed in Fig S24E UMAP figure panel using Illustrator
pdf('/filepath/fig5/fig5_supp/flu_tfh_gng4_feat_plot_with_legend.pdf', width = 6, height = 6)
FeaturePlot(flu_tfh_obj, features = 'GNG4', pt.size = 2, order = TRUE, cols = c('grey','coral2')) + coord_fixed() + theme(plot.title = element_blank())
dev.off()

# Fig S24F - DotPlot of RNA expression for GNG4 and other DEGs across clusters of L2 Tfh subclustering object (reanalysis of Schattgen et al. study) ----

# Set cluster order for dotplot
tfh_ident_order <- c('cycling','naive','pre/memory','Treg','IL10 TFH','GC')
Idents(flu_tfh_obj) <- 'Tfh_type'
Idents(flu_tfh_obj) <- factor(
  Idents(flu_tfh_obj),
  levels = tfh_ident_order
)

# Make dotplot
pdf('/filepath/fig5/fig5_supp/flu_tfh_feat_dotplot.pdf', height = 5, width = 10)
DotPlot(flu_tfh_obj, features = c('MKI67','LEF1','CCR7','IL7R','KLF2','GPR183','FOXP3','IL2RA','IL10','CXCR5','PDCD1','TIGIT','BCL6','TOX2','GNG4'), dot.scale = 10) + scale_color_gradient2(
  low = rdbu_colors[11],
  mid = "white",
  high = rdbu_colors[2],
  midpoint = 0
) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'italic', family = 'sans', angle = 45, hjust = 1))
dev.off()

# Fig S24G - UMAP of Level 1 T cell object from Schattgen et al. study ----

Idents(flu_tcell_obj) <- 'seurat_clusters' # original author annotations, where c14 is the original 'Tfh/Treg' cluster (no 'GC Tfh' annotation as is)
pdf('/filepath/fig5/fig5_supp/flu_tcell_clust_numb_dimplot.pdf', width = 6, height = 6)
DimPlot(flu_tcell_obj, group.by = 'seurat_clusters', cols = cols_schattgen, order = TRUE, pt.size = 2, label.box = TRUE, label = TRUE) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()

# Fig S24H -  UMAP of GNG4 RNA expression in L1 T cell object (reanalysis of Schattgen et al. study)  ----

# UMAP without legend
pdf('/filepath/fig5/fig5_supp/flu_tcell_no_labs_gng4_feat_plot.pdf', width = 6, height = 6)
FeaturePlot(flu_tcell_obj, features = 'GNG4', pt.size = 2, label = FALSE, order = TRUE, cols = c('grey','coral2')) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()

# UMAP with legend, then cropped to legend and placed in Fig S24H UMAP figure panel using Illustrator
Idents(flu_tcell_obj) <- 'seurat_clusters' # where c14 is the original 'Tfh/Treg' cluster
pdf('/filepath/fig5/fig5_supp/flu_tcell_gng4_feat_plot_legend.pdf', width = 6, height = 6)
FeaturePlot(flu_tcell_obj, features = 'GNG4', pt.size = 2, label = TRUE, order = TRUE, cols = c('grey','coral2')) + coord_fixed() + theme(plot.title = element_blank())
dev.off()

# Fig S24I - DotPlot of RNA expression for GNG4 and other DEGs across clusters of L1 T cell object (reanalysis of Schattgen et al. study) ----

pdf('/filepath/fig5/fig5_supp/flu_tcell_seurat_cluster_dotplot.pdf', height = 5, width = 7.75)
Idents(flu_tcell_obj) <- 'seurat_clusters' # where c14 is the original 'Tfh/Treg' cluster
DotPlot(flu_tcell_obj, features = c('CD4','CD8A','MKI67','LEF1','CCR7','IL7R','KLF2','GPR183','FOXP3','IL2RA','IL10','CXCR5','PDCD1','TIGIT','BCL6','TOX2','GNG4'), dot.scale = 6,
        group.by = 'seurat_clusters', cluster.idents = TRUE) +
  scale_color_gradient2(
    low = rdbu_colors[11],
    mid = "white",
    high = rdbu_colors[2],
    midpoint = 0          
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'italic', family = 'sans', angle = 45, hjust = 1))
dev.off()

# Fig S24J - UMAP of L2-to-L1 mapped 'GC Tfh' cluster versus other clusters of L1 T cell object (reanalysis of Schattgen et al. study)  -----

Idents(flu_tcell_obj) <- 'gctfh_annot' # from metadata column determined above for mapping L2 cell annotations to L1 object
gctfh_yesno_cols <- c('no' = 'grey', # no, not from 'GC' cluster in L2 object
                      'yes' = '#e17794') # yes, is from 'GC' cluster in L2 object
pdf('/filepath/fig5/fig5_supp/flu_gctfh_mapped_to_tcell_obj_annot_dimplot.pdf', width = 6, height = 6)
DimPlot(flu_tcell_obj, group.by = 'gctfh_annot', cols = gctfh_yesno_cols, order = TRUE, pt.size = 2) + coord_fixed() + theme(plot.title = element_blank()) + NoLegend()
dev.off()

# Fig S24K - DotPlot of RNA expression for GNG4 and other DEGs between L2-to-L1 mapped 'GC Tfh' cluster versus other clusters of L1 T cell object (reanalysis of Schattgen et al. study)  -----
Idents(flu_tcell_obj) <- 'gctfh_vs_others_annot' # metadata column determined above wherein L2 'GC Tfh' cells are annotated in L1 object and original L1 'seurat_clusters' for remaining cells in object are retained
pdf('/filepath/fig5/fig5_supp/flu_tcell_seurat_clust_vs_gctfh_dotplot.pdf', height = 5, width = 8.25)
DotPlot(flu_tcell_obj, features = c('CD4','CD8A','MKI67','LEF1','CCR7','IL7R','KLF2','GPR183','FOXP3','IL2RA','IL10','CXCR5','PDCD1','TIGIT','BCL6','TOX2','GNG4'), dot.scale = 5.25,
        group.by = 'gctfh_vs_others_annot', cluster.idents = TRUE) +
  scale_color_gradient2(
    low = rdbu_colors[11],
    mid = "white",
    high = rdbu_colors[2],
    midpoint = 0
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'italic', family = 'sans', angle = 45, hjust = 1))
dev.off()