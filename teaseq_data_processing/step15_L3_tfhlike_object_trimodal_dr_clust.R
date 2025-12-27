# Step 15 - Level 3 TEAseq Analysis - Subclustering of Tfh-like clusters from L2 T cell object, as annotated in Step 14

# 1) Set up R working environment ----

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

# Set working directory for this step
setwd("/filepath/step15_tfh_subcluster/")

# Set color palette
cols_numb <- c(
  '0' = "#e49a78", '1' = "#d66d97", '2' = "#80588b", '3' = "#90a6bf", '4' = "#9ce4cc",
  '5' = "#4b94d6", '6' = "#f4a6c7", '7' = "#d4b5e4", '8' = "#6e485a", '9' = "#50664f",
  '10' = "#7b7f9c", '11' = "#723d56", '12' = "#b3f5dd", '13' = "#8ab1cc", '14' = "#535c81",
  '15' = "#e1d4c7", '16' = "#c9a3c1", '17' = "#add5df", '18' = "#e17794", '19' = "#72d1b4",
  '20' = "#94a890", '21' = "#dff3ed", '22' = "#4a3c5d", '23' = "#e0c1b4", '24' = "#f7f1e0"
)

# 2) Import Seurat object from TEAseq Data Preprocessing Step 14 (trimodal dimensionality reduction and L2 3WNN subclustering of T cells from L1 object, including Harmony integration across donors) ----
l2_teaseq_tcell_obj <- readRDS('/step14_l2_teaseq_tcell_analysis/l2_teaseq_tcell_obj.rds')

# 3) Filter T cell object for CXCR5+ Tfh-like clusters ----

# Three clusters (c4, c5, and c17) were found to be CXCR5-enriched Tfh-like clusters in Step 14 T Cell subclustering step
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'
l3_teaseq_tfh_obj <- subset(l2_teaseq_tcell_obj, idents = c('Tfh GC', 'Tfh IL10', 'CD4 Tcm/fh'))

# 4) ATAC - dimensionality reduction before integration ----

# Standard Signac workflow - excluding first dim due to strong correlation with ATAC read depth and tuned to dims 2:10
DefaultAssay(l3_teaseq_tfh_obj) <- "ATAC"
l3_teaseq_tfh_obj <- RunTFIDF(l3_teaseq_tfh_obj, assay = 'ATAC')
l3_teaseq_tfh_obj <- FindTopFeatures(l3_teaseq_tfh_obj, min.cutoff = 'q0', assay = 'ATAC') # using all features
l3_teaseq_tfh_obj <- RunSVD(l3_teaseq_tfh_obj, reduction.name = "lsi.atac", reduction.key = "atacLSI_", assay = 'ATAC')

# Evaluating correlation of dims with ATAC read depth
DepthCor(l3_teaseq_tfh_obj, reduction = 'lsi.atac', assay = 'ATAC')
ElbowPlot(l3_teaseq_tfh_obj, ndims = 50, reduction = "lsi.atac")
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'atacLSI_1', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'atacLSI_2', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'atacLSI_3', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'atacLSI_4', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'atacLSI_5', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))

# Standard Seurat workflow
l3_teaseq_tfh_obj <- FindNeighbors(l3_teaseq_tfh_obj, reduction= "lsi.atac", dims = 2:11, assay = 'ATAC')
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, reduction = 'lsi.atac', dims = 2:11, reduction.name = "umap.atac", reduction.key = "atacUMAP_", assay = 'ATAC')
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, resolution = 0.4, cluster.name = "atac.tfh.subcluster", graph.name = "ATAC_snn")

# Visualization
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('ATAC - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('ATAC - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('ATAC - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('ATAC - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('ATAC - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('ATAC - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('ATAC - Before Harmony') + coord_fixed()

# Feature exploration
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster'
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, features = 'rna_CXCR5', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()

# 5) ATAC - Harmony integration ----

# ATAC data processed as a Signac 'ChromatinAssay' - unlike RNA/ATAC no layers to split or join
l3_teaseq_tfh_obj <- RunHarmony(
  object = l3_teaseq_tfh_obj,
  group.by.vars = 'hto.donor', # integrating across donors
  assay.use = 'ATAC',
  reduction.use = 'lsi.atac',
  reduction.save = 'lsi.atac.harmony',
  project.dim = FALSE,
  dims.use = 2:11,
  seed = 26
)

# Checking whether first LSI is still primarily explained by ATAC read depth
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'lsiatacharmony_1', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'lsiatacharmony_2', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'lsiatacharmony_3', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'lsiatacharmony_4', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
FeatureScatter(l3_teaseq_tfh_obj, feature1 = 'lsiatacharmony_5', feature2 = 'nCount_ATAC', group.by = 'hto.tissue', cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6'))
ElbowPlot(l3_teaseq_tfh_obj, reduction = 'lsi.atac.harmony', ndims = 10)

# Standard Seurat workflow with resulting 10 dims after Harmony integration of dims 2:11
l3_teaseq_tfh_obj <- FindNeighbors(l3_teaseq_tfh_obj, reduction = "lsi.atac.harmony", dims = 1:10, assay = 'ATAC', graph.name = c('ATAC_harmony_nn','ATAC_harmony_snn'))
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, dims = 1:10, reduction = "lsi.atac.harmony", reduction.name = "umap.atac.harmony", reduction.key = "atacharmonyUMAP_", assay = 'ATAC')
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, resolution = 0.4, cluster.name = "atac.tfh.subcluster.harmony", assay = 'ATAC', graph.name = 'ATAC_harmony_snn')

# Visualization
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('ATAC - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('ATAC - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('ATAC - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('ATAC - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('ATAC - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('ATAC - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('ATAC - After Harmony') + coord_fixed()

# Feature exploration
Idents(l3_teaseq_tfh_obj) <- 'atac.tfh.subcluster.harmony'
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'nCount_ATAC', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()

# 6) RNA - dimensionality reduction before integration ----

# Standard Seurat workflow
DefaultAssay(l3_teaseq_tfh_obj) <- "RNA"
l3_teaseq_tfh_obj <- NormalizeData(l3_teaseq_tfh_obj, assay = 'RNA')
l3_teaseq_tfh_obj <- FindVariableFeatures(l3_teaseq_tfh_obj, assay = 'RNA')
l3_teaseq_tfh_obj <- ScaleData(l3_teaseq_tfh_obj, assay = 'RNA')
l3_teaseq_tfh_obj <- RunPCA(l3_teaseq_tfh_obj, assay = 'RNA', reduction.name = 'pca.rna', reduction.key = 'rnaPC_')
ElbowPlot(l3_teaseq_tfh_obj, ndims = 50, reduction = 'pca.rna') + ggtitle('RNA Elbow Plot - Before Harmony')
l3_teaseq_tfh_obj <- FindNeighbors(l3_teaseq_tfh_obj, assay = 'RNA', reduction = 'pca.rna', dims = 1:35)
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, assay = 'RNA', reduction = 'pca.rna', dims = 1:35, reduction.name = 'umap.rna')
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, assay = 'RNA', cluster.name = 'rna.tfh.subcluster', resolution = 0.7, graph.name = 'RNA_snn')

# Clustering visualization before integration
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.sort'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, group.by = 'hto.sort', pt.size = 1.5) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, split.by = 'hto.sort', pt.size = 1.5, cols = cols_numb, ncol = 4) + ggtitle('RNA - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, split.by = 'hto.donor', pt.size = 1.5, cols = cols_numb, ncol = 4) + ggtitle('RNA - Before Harmony') + coord_fixed()

# Feature exploration
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('lightgrey', 'red4')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'rna_CD4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'rna_CD8A', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'rna_CXCR5', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()

# 7) RNA - Harmony integration ----

# Integrating RNA layers
l3_teaseq_tfh_obj[["RNA"]] <- split(l3_teaseq_tfh_obj[["RNA"]], f = l3_teaseq_tfh_obj$hto.donor) # split RNA layers by donor
l3_teaseq_tfh_obj <- IntegrateLayers(object = l3_teaseq_tfh_obj, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca.rna", 
                                 new.reduction = "pca.rna.harmony",
                                 assay = 'RNA',
                                 verbose = TRUE,
                                 seed = 26)

# Standard Seurat workflow - tuned to 30 RNA Harmony dims
l3_teaseq_tfh_obj <- FindNeighbors(l3_teaseq_tfh_obj, reduction = "pca.rna.harmony", dims = 1:50, graph.name = c('RNA_harmony_nn','RNA_harmony_snn'))
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, dims = 1:50, reduction = "pca.rna.harmony", reduction.name = "umap.rna.harmony", reduction.key = "rnaharmonyUMAP_")
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, resolution = 0.7, cluster.name = "rna.tfh.subcluster.harmony", graph.name = 'RNA_harmony_snn')

# Visualization after integration
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('RNA - After Harmony')
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('RNA - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('RNA - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('RNA - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('RNA - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('RNA - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('RNA - After Harmony') + coord_fixed()

# Feature exploration
Idents(l3_teaseq_tfh_obj) <- 'rna.tfh.subcluster'
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_CD4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_BCL6', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_CXCR5', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()


# 8) ADT - dimensionality reduction before integration ----

# Standard Seurat workflow - tuned to 25 PC for ADT
DefaultAssay(l3_teaseq_tfh_obj) <- "ADT"
l3_teaseq_tfh_obj <- NormalizeData(l3_teaseq_tfh_obj, assay= "ADT", normalization.method = "CLR", margin = 2)
l3_teaseq_tfh_obj <- ScaleData(l3_teaseq_tfh_obj, assay = "ADT")
l3_teaseq_tfh_obj <- FindVariableFeatures(l3_teaseq_tfh_obj, assay = "ADT")
l3_teaseq_tfh_obj <- RunPCA(l3_teaseq_tfh_obj, assay = "ADT", reduction.name = "pca.adt", reduction.key = "adtPC_")
ElbowPlot(l3_teaseq_tfh_obj, ndims = 50, reduction = "pca.adt")
l3_teaseq_tfh_obj <- FindNeighbors(l3_teaseq_tfh_obj, reduction = "pca.adt", dims = 1:25)
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, assay = "ADT", dims = 1:25, reduction = "pca.adt", reduction.name = "umap.adt", reduction.key = "adtUMAP_")
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, resolution = 1, cluster.name = "adt.tfh.subcluster", graph.name = "ADT_snn")

# Visualization
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('ADT - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('ADT - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('ADT - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('ADT - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('ADT - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('ADT - Before Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('ADT - Before Harmony') + coord_fixed()

# Feature exploration
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster'
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'rna_CD4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'rna_CD8A', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'rna_CXCR5', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()

# 9) ADT - Harmony integration ----

# Integrating ADT layers across donors
l3_teaseq_tfh_obj[["ADT"]] <- split(l3_teaseq_tfh_obj[["ADT"]], f = l3_teaseq_tfh_obj$hto.donor) # split ADT layers by donor
l3_teaseq_tfh_obj <- IntegrateLayers(object = l3_teaseq_tfh_obj, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca.adt", 
                                 new.reduction = "pca.adt.harmony",
                                 assay = 'ADT',
                                 verbose = TRUE,
                                 seed = 26)

# Standard Seurat workflow - tuned to 25 Harmony ADT dims
l3_teaseq_tfh_obj <- FindNeighbors(l3_teaseq_tfh_obj, reduction = "pca.adt.harmony", dims = 1:25, assay = 'ADT', graph.name = c('ADT_harmony_nn','ADT_harmony_snn'))
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, dims = 1:25, reduction = "pca.adt.harmony", reduction.name = "umap.adt.harmony", reduction.key = "adtharmonyUMAP_", assay = 'ADT')
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, resolution = 0.8, cluster.name = "adt.tfh.subcluster.harmony", assay = 'ADT', graph.name = 'ADT_harmony_snn')

# Visualization after integration
Idents(l3_teaseq_tfh_obj) <- "adt.tfh.subcluster.harmony"
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('ADT - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('ADT - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('ADT - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('ADT - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('ADT - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('ADT - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('ADT - After Harmony') + coord_fixed()

# Feature exploration
Idents(l3_teaseq_tfh_obj) <- 'adt.tfh.subcluster.harmony'
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_IL10', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_IL7R', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'nCount_ATAC', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()

# 10) 3WNN clustering and dimensionality reduction ----

# Using same number of Harmony dims as tuned above for each assay - 35 RNA, 10 ATAC (following integration of LSI components 2:11), and 25 ADT
l3_teaseq_tfh_obj <- FindMultiModalNeighbors(l3_teaseq_tfh_obj, 
                                       reduction.list = list("pca.rna.harmony", "lsi.atac.harmony", "pca.adt.harmony"), 
                                       dims.list = list(1:35, 1:10, 1:25),
                                       knn.graph.name = 'knn.wnn.harmony',
                                       snn.graph.name = 'snn.wnn.harmony',
                                       weighted.nn.name = 'wnn.harmony'
                                       )
l3_teaseq_tfh_obj <- RunUMAP(l3_teaseq_tfh_obj, nn.name = "wnn.harmony", reduction.name = "umap.wnn.harmony", reduction.key = "harmonywnnUMAP_")
l3_teaseq_tfh_obj <- FindClusters(l3_teaseq_tfh_obj, graph.name = "snn.wnn.harmony", algorithm = 3, resolution = 0.6, verbose = TRUE, cluster.name = "wnn.tfh.subcluster.harmony")

# Visualization
Idents(l3_teaseq_tfh_obj) <- "wnn.tfh.subcluster.harmony"
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = TRUE, cols = cols_numb, pt.size = 1.5, label.box = TRUE) + ggtitle('3WNN - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'asab'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, group.by = 'asab', pt.size = 1.5, cols = c('Male' = '#4b94d6', 'Female' = '#e49a78')) + ggtitle('3WNN - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'wnn.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, split.by = 'asab', pt.size = 1.5, cols = cols_numb) + ggtitle('3WNN - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'hto.tissue'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, group.by = 'hto.tissue', pt.size = 1.5, cols = c('Tonsil' = '#e49a78', 'PBMC' = '#4b94d6')) + ggtitle('3WNN - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'wnn.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, split.by = 'hto.tissue', pt.size = 1.5, cols = cols_numb) + ggtitle('3WNN - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'CC.Phase'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, group.by = 'CC.Phase', cols = c('S' = '#4b94d6', 'G1' = '#e49a78', 'G2M' = '#723d56'), pt.size = 1.5) + ggtitle('3WNN - After Harmony') + coord_fixed()
Idents(l3_teaseq_tfh_obj) <- 'wnn.tfh.subcluster.harmony'
DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, split.by = 'CC.Phase', cols = cols_numb, pt.size = 1.5) + ggtitle('3WNN - After Harmony') + coord_fixed()

# Feature exploration
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'rna_GNG4', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'rna_IL7R', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'percent.ribo', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'percent.mt', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()
FeaturePlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'nCount_ATAC', pt.size = 1.5, order = TRUE, cols = c('#D3D3D3','#7F1E2B')) + coord_fixed()

# 11) 3WNN subcluster annotation ----

# Assign names to Tfh subcluster numbers
tfh_subcluster_names <- c(
  "0" = "Tfh-Circ",
  "1" = "Tcm",
  "2" = "Tfh-Int",
  "3" = "Tfh-NFATC1",
  "4" = "Tfh-IL10",
  "5" = "Tfh-Resting",
  "6" = "Tfh-CXCL13",
  "7" = "Tfh-BOB1",
  "8" = "Tfh-AP1"
)

# Rename idents
l3_teaseq_tfh_obj$tfh_wnn_annot <- l3_teaseq_tfh_obj$wnn.tfh.subcluster.harmony
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
l3_teaseq_tfh_obj <- RenameIdents(l3_teaseq_tfh_obj, tfh_subcluster_names)
l3_teaseq_tfh_obj$tfh_wnn_annot <- Idents(l3_teaseq_tfh_obj)
table(Idents(l3_teaseq_tfh_obj))


# Add GC vs nonGC-like group annotation
gc_vs_nongc_like <- c(
  "Tfh-Circ" = 'nonGC',
  "Tcm" = "nonGC",
  "Tfh-Int" = "GC",
  "Tfh-NFATC1" = "GC",
  "Tfh-IL10" = "GC",
  "Tfh-Resting" = "nonGC",
  "Tfh-CXCL13" = "GC",
  "Tfh-BOB1" = "GC",
  "Tfh-AP1" = "nonGC"
)
l3_teaseq_tfh_obj$gc_vs_nongc_like <- l3_teaseq_tfh_obj$tfh_wnn_annot
Idents(l3_teaseq_tfh_obj) <- 'gc_vs_nongc_like'
l3_teaseq_tfh_obj <- RenameIdents(l3_teaseq_tfh_obj, gc_vs_nongc_like)
l3_teaseq_tfh_obj$gc_vs_nongc_like <- Idents(l3_teaseq_tfh_obj)
table(l3_teaseq_tfh_obj$gc_vs_nongc_like)

# Visualization
tfh.colors.assigned <- c('Tfh-Circ' = "#e49a78", 'Tcm' = "#d66d97", 'Tfh-Int' = "#80588b", 'Tfh-NFATC1' = "#90a6bf", 'Tfh-IL10' = "#9ce4cc",'Tfh-Resting' = "#4b94d6", 'Tfh-CXCL13' = "#f4a6c7", 'Tfh-BOB1' = "#d4b5e4", 'Tfh-AP1' = "#6e485a")
do_DimPlot(l3_teaseq_tfh_obj, reduction = "umap.wnn.harmony", group.by = 'tfh_wnn_annot', label = TRUE, label.size = 4, pt.size = 2, repel = FALSE, font.size = 9, colors.use = tfh.colors.assigned, border.size = 1.25) + coord_fixed() + NoLegend()

# 12) Add TF motif information to object ----

# Guided by Signac vignette - https://stuartlab.org/signac/articles/motif_vignette

# Filtering feature names - https://github.com/stuart-lab/signac/issues/780
DefaultAssay(l3_teaseq_tfh_obj) <- 'ATAC'
gr <- granges(l3_teaseq_tfh_obj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
l3_teaseq_tfh_obj[['ATAC']] <- subset(l3_teaseq_tfh_obj[["ATAC"]], features = feat.keep)

# Get motif PFM information
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

l3_teaseq_tfh_obj <- AddMotifs(
  object = l3_teaseq_tfh_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Add ATAC ChromVAR assay to object ----
l3_teaseq_tfh_obj <- RunChromVAR(
  object = l3_teaseq_tfh_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = 'ATAC',
  new.assay.name = 'chromvar'
)

# 13) Add ATAC GeneActivity assay information to object ----

# Modify parameters to include all gene biotypes and widths - will enable closer comparison of ATAC and larger RNA feature sets
DefaultAssay(l3_teaseq_tfh_obj) <- 'ATAC'
l3_teaseq_geneactivity_mtx <- GeneActivity(
  l3_teaseq_tfh_obj,
  assay = 'ATAC',
  features = rownames(l3_teaseq_tfh_obj@assays$RNA), # consider all genes in RNA assay
  extend.upstream = 2000, # default
  extend.downstream = 0, # default
  biotypes = NULL, # include noncoding genes, otherwise default 'protein_coding'
  max.width = NULL, # do not filter based on gene length, otherwise default '5e+05'
  process_n = 2000, # decreasing runtime, default 2000
  gene.id = FALSE,
  verbose = TRUE
)

# Create GeneActivity assay
l3_teaseq_tfh_obj[['ACT']] <- CreateAssayObject(counts = l3_teaseq_geneactivity_mtx)

# Normalize following Signac tutorial
l3_teaseq_tfh_obj <- NormalizeData(
  object = l3_teaseq_tfh_obj,
  assay = 'ACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(l3_teaseq_tfh_obj$nCount_ACT)
)

# 14) Add SCENIC assay to object ----

setwd('/filepath/l3_tfh_scenic_output_folder/')
scenicOptions <- readRDS("/filepath/l3_tfh_scenic_output_folder/scenicOptions4.rds") # after complete pipeline run, refer to SCENIC code folder
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC") # load AUC (regulon activity) matrix
AUCmatrix <- getAUC(regulonAUC)
l3_teaseq_tfh_obj[["SCENIC"]] <- CreateAssayObject(data = AUCmatrix)

# 15) Find differentially expressed features between L3 Tfh-like subclusters across all modalities----

# Set cluster identity to 3WNN L3 cluster names
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'

# Join RNA and ADT layers to run FindAllMarkers - other assays do not have layers split by donor
l3_teaseq_tfh_obj <- JoinLayers(l3_teaseq_tfh_obj, assay = 'RNA')
l3_teaseq_tfh_obj <- JoinLayers(l3_teaseq_tfh_obj, assay = 'ADT')

# Find TEAseq 3WNN L3 Tfh-like cluster markers for each assay using FindAllMarkers, adjusted as detailed below to fit the data type of each assay
# ATAC - peaks ('ATAC'), GeneActivity ('ACT', promoter region and gene body accessibility), and chromVAR (TF motif deviation Z-scores)
# RNA - transcripts ('RNA', log-normalized) and SCENIC (TF regulon AUC scores)
# ADT - antibody-derived tag epitope information (CLR-normalized)

# TEAseq 3WNN L3 Tfh-like Cluster Markers - ATAC (peaks)
DefaultAssay(l3_teaseq_tfh_obj) <- 'ATAC'
l3_tfh_wnn_clust_atac_peak_markers <- FindAllMarkers(l3_teaseq_tfh_obj, assay = 'ATAC') # using default parameters
dap_names <- l3_tfh_wnn_clust_atac_peak_markers$gene
dap_gr <- StringToGRanges(dap_names)
dap_closest_gene <- ClosestFeature(
  object = l3_teaseq_tfh_obj,
  regions    = dap_gr,
  annotation = Annotation(l3_teaseq_tfh_obj)
)
l3_tfh_wnn_clust_atac_peak_markers <- l3_tfh_wnn_clust_atac_peak_markers %>%
  mutate(
    closest_gene_symbol = dap_closest_gene$gene_name,
    closest_ensembl_id = dap_closest_gene$gene_id
  )
l3_tfh_wnn_clust_atac_peak_markers <- l3_tfh_wnn_clust_atac_peak_markers %>%
  rename(
    l3_tfh_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    atac_peak = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l3_tfh_wnn_clust_atac_peak_markers <- l3_tfh_wnn_clust_atac_peak_markers %>% arrange(l3_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l3_tfh_wnn_clust_atac_peak_markers <- l3_tfh_wnn_clust_atac_peak_markers %>% relocate(c('l3_tfh_cluster','atac_peak','closest_gene_symbol','closest_ensembl_id','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l3_tfh_wnn_clust_atac_peak_markers, '/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_atac_peak_markers.rds')

# TEAseq 3WNN L3 Tfh-like Cluster Markers - ACT (ATAC-based Signac GeneActivity summary of gene body and promoter region accessibility)
DefaultAssay(l3_teaseq_tfh_obj) <- 'ACT'
l3_tfh_wnn_clust_atac_geneactivity_markers <- FindAllMarkers(l3_teaseq_tfh_obj, assay = 'ACT') # using default parameters
l3_tfh_wnn_clust_atac_geneactivity_markers <- l3_tfh_wnn_clust_atac_geneactivity_markers %>%
  rename(
    l3_tfh_cluster = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l3_tfh_wnn_clust_atac_geneactivity_markers <- l3_tfh_wnn_clust_atac_geneactivity_markers %>% arrange(l3_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l3_tfh_wnn_clust_atac_geneactivity_markers <- l3_tfh_wnn_clust_atac_geneactivity_markers %>% relocate(c('l3_tfh_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l3_tfh_wnn_clust_atac_geneactivity_markers, '/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_atac_geneactivity_markers.rds')

# TEAseq 3WNN L3 Tfh-like Cluster Markers - chromVAR (ATAC-based TF motif deviation Z-scores)
l3_tfh_wnn_clust_atac_chromvar_markers <- FindAllMarkers(l3_teaseq_tfh_obj, assay = 'chromvar', 
                                                         mean.fxn = rowMeans, # use row means function to find cluster differences as chromVAR Z-scores are not log-normalized like Seurat RNA data
                                                         fc.name = "mean_diff_z_score") # using default parameters
pfm <- getMatrixSet( # get TF motif PFM information
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
id_to_tf <- sapply(pfm, function(x) name(x)) # function to convert PFM ID names to common TF names
names(id_to_tf) <- sapply(pfm, ID)
l3_tfh_wnn_clust_atac_chromvar_markers$tf_name <- id_to_tf[l3_tfh_wnn_clust_atac_chromvar_markers$gene]
l3_tfh_wnn_clust_atac_chromvar_markers <- l3_tfh_wnn_clust_atac_chromvar_markers %>%
  rename(
    l3_tfh_cluster = cluster,
    jaspar_pfm_id = gene,
    p_val_raw = p_val,
    pct_nonzero_deviation_in_clust = pct.1, # rename, as nonzero chromVAR Z-scores do not mean 'positive' as with RNA expression 
    pct_nonzero_deviation_in_others = pct.2
  )
l3_tfh_wnn_clust_atac_chromvar_markers <- l3_tfh_wnn_clust_atac_chromvar_markers %>% arrange(l3_tfh_cluster, desc(mean_diff_z_score), p_val_adj)
l3_tfh_wnn_clust_atac_chromvar_markers <- l3_tfh_wnn_clust_atac_chromvar_markers %>% relocate(c('l3_tfh_cluster','tf_name','jaspar_pfm_id','mean_diff_z_score','p_val_raw','p_val_adj'))
saveRDS(l3_tfh_wnn_clust_atac_chromvar_markers, '/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_atac_chromvar_markers.rds')

# TEAseq 3WNN L3 Tfh-like Cluster Markers - RNA (standard Seurat log-normalized mRNA count data)
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
l3_tfh_wnn_clust_rna_markers <- FindAllMarkers(l3_teaseq_tfh_obj, assay = 'RNA') # using default parameters
l3_tfh_wnn_clust_rna_markers <- l3_tfh_wnn_clust_rna_markers %>%
  rename(
    l3_tfh_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l3_tfh_wnn_clust_rna_markers <- l3_tfh_wnn_clust_rna_markers %>% arrange(l3_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l3_tfh_wnn_clust_rna_markers <- l3_tfh_wnn_clust_rna_markers %>% relocate(c('l3_tfh_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l3_tfh_wnn_clust_rna_markers, '/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_rna_markers.rds')

# TEAseq 3WNN L3 Tfh-like Cluster Markers - SCENIC (RNA-based TF regulon inference)
DefaultAssay(l3_teaseq_tfh_obj) <- "SCENIC"
l3_tfh_wnn_clust_rna_scenic_markers <- FindAllMarkers(l3_teaseq_tfh_obj, assay = 'SCENIC', 
                                                      mean.fxn = rowMeans, # SCENIC regulon AUC scores are not log-transformed
                                                      fc.name = "mean_diff_AUC",
                                                      min.pct = 0, logfc.threshold = 0) # adjust default filters as the mean difference in AUC may be lower than 0.1 
l3_tfh_wnn_clust_rna_scenic_markers <- l3_tfh_wnn_clust_rna_scenic_markers %>%
  rename(
    l3_tfh_cluster = cluster,
    scenic_regulon = gene,
    p_val_raw = p_val,
    pct_nonzero_auc_in_clust = pct.1,
    pct_nonzero_auc_others = pct.2
  )
l3_tfh_wnn_clust_rna_scenic_markers <- l3_tfh_wnn_clust_rna_scenic_markers %>% arrange(l3_tfh_cluster, desc(mean_diff_AUC), p_val_adj)
l3_tfh_wnn_clust_rna_scenic_markers <- l3_tfh_wnn_clust_rna_scenic_markers %>% relocate(c('l3_tfh_cluster','scenic_regulon','mean_diff_AUC','p_val_raw','p_val_adj'))
saveRDS(l3_tfh_wnn_clust_rna_scenic_markers, '/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_rna_scenic_markers.rds')

# TEAseq 3WNN L3 Tfh-like Cluster Markers - ADT (CLR-normalized count data)
DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
l3_tfh_wnn_clust_adt_markers <- FindAllMarkers(l3_teaseq_tfh_obj, assay = 'ADT') # using default parameters
l3_tfh_wnn_clust_adt_markers <- l3_tfh_wnn_clust_adt_markers %>%
  rename(
    l3_tfh_cluster = cluster,
    adt_epitope = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l3_tfh_wnn_clust_adt_markers <- l3_tfh_wnn_clust_adt_markers %>% arrange(l3_tfh_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l3_tfh_wnn_clust_adt_markers$adt_epitope <- gsub("^adt-", "", l3_tfh_wnn_clust_adt_markers$adt_epitope)
l3_tfh_wnn_clust_adt_markers <- l3_tfh_wnn_clust_adt_markers %>% relocate(c('l3_tfh_cluster','adt_epitope','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l3_tfh_wnn_clust_adt_markers, '/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_adt_markers.rds')

# Export all cluster marker files as sheets within one XLSX spreadsheet
clust_marker_df_dir   <- '/filepath/step15_tfh_subcluster/teaseq_l3_tfh_cluster_markers'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) |> str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/step15_tfh_subcluster/l3_tfh_wnn_clust_markers_all_assays.xlsx")  

# 16) Save L3 object and proceed to main analysis and data visualization steps ----
saveRDS(l3_teaseq_tfh_obj, '/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')