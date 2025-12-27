# Step 14 - Level 2 TEAseq Analysis - Subclustering of T cell clusters from L1 object, as annotated in Step 13
# Trimodal dimensionality reduction, Harmony integration, weighted-nearest neighbor analysis, clustering, and annotation

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

# Set working directory
setwd("/filepath/step14_l2_teaseq_tcell_analysis/")

# Create color palette for visualization
tcell_cols_unspec <- c(
  "#8ab1cc", "#add5df", "#94a890", "#e0c1b4", "#e49a78",
  "#657ab0", "#e17794", "#72d1b4", "#c9744d", "#d66d97",
  "#50664f", "#d4b5e4", "#e1d4c7", "#c9a3c1", "#f4a6c7",
  "#bad4f4", "#dfcc78", "#728762", "#c49980", "#b6b7ba",
  "#d699ac", "#50664f", "#4a3c5d", "#69c9d1", "#f7f1e0",
  "#c9744d", "#b37cbf")

# Numbers assigned colors
tcell_cols_numb <- c(
  "0"  = "#8ab1cc",
  "1"  = "#add5df",
  "2"  = "#94a890",
  "3"  = "#e0c1b4",
  "4"  = "#e49a78",
  "5"  = "#657ab0",
  "6"  = "#e17794",
  "7"  = "#72d1b4",
  "8"  = "#c9744d",
  "9"  = "#d66d97",
  "10" = "#50664f",
  "11" = "#d4b5e4",
  "12" = "#e1d4c7",
  "13" = "#c9a3c1",
  "14" = "#f4a6c7",
  "15" = "#bad4f4",
  "16" = "#dfcc78",
  "17" = "#728782",
  "18" = "#c49980",
  "19" = "#b6b7ba",
  "20" = "#d699ac",
  "21" = "#728762",
  "22" = "#4a3c5d",
  "23" = "#69c9d1",
  "24" = "#f7f1e0",
  "25" = "#c9744d",
  "26" = "#b37cbf"
)

# Tissue colors
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# 2) Import Seurat object from TEAseq Data Preprocessing Step 13 (trimodal dimensionality reduction and L1 3WNN clustering of all tonsil and peripheral blood mononuclear cells, including Harmony integration across donors) ----
l1_teaseq_obj <- readRDS('/step13_l1_teaseq_analysis/l1_teaseq_obj.rds')

# 3) Filter bulk object for annotated T Cell clusters ----
l2_teaseq_tcell_obj <- subset(l1_teaseq_obj, idents = c('CD4 Tn','CD4 Tm','Treg','Tfh-like','CD8 Tn','CD8 Tm/Inn','CD8 UTC')) # Excluding L1 'CD4 Tn Ribo' cluster, as noted in supplementary methods
ncol(l2_teaseq_tcell_obj)

# 4) ATAC - dimensionality reduction before integration ----

# Standard ATAC DR and clustering workflow from Signac
DefaultAssay(l2_teaseq_tcell_obj) <- "ATAC"
l2_teaseq_tcell_obj <- RunTFIDF(l2_teaseq_tcell_obj, assay = 'ATAC')
l2_teaseq_tcell_obj <- FindTopFeatures(l2_teaseq_tcell_obj, min.cutoff = 'q0', assay = 'ATAC')
l2_teaseq_tcell_obj <- RunSVD(l2_teaseq_tcell_obj, reduction.name = "lsi.atac", reduction.key = "atacLSI_", assay = 'ATAC')

# Inspecting LSI components for correlation with ATAC read depth
ElbowPlot(l2_teaseq_tcell_obj, ndims = 50, reduction = "lsi.atac")
FeatureScatter(l2_teaseq_tcell_obj, feature1 = 'atacLSI_1', feature2 = 'nCount_ATAC')
FeatureScatter(l2_teaseq_tcell_obj, feature1 = 'atacLSI_2', feature2 = 'nCount_ATAC')
FeatureScatter(l2_teaseq_tcell_obj, feature1 = 'atacLSI_3', feature2 = 'nCount_ATAC')
FeatureScatter(l2_teaseq_tcell_obj, feature1 = 'atacLSI_4', feature2 = 'nCount_ATAC')

# Excluding first LSI component per Signac tutorial guidelines
l2_teaseq_tcell_obj <- FindNeighbors(l2_teaseq_tcell_obj, reduction= "lsi.atac", dims = 2:50, assay = 'ATAC')
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, reduction = 'lsi.atac', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_", assay = 'ATAC')
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, resolution = 0.5, cluster.name = "atac.tcell.subcluster", graph.name = "ATAC_snn")
Idents(l2_teaseq_tcell_obj) <- 'atac.tcell.subcluster'

# Clustering visualization before harmony integration
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.atac", label = TRUE) + ggtitle('ATAC - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.atac", label = TRUE, split.by = 'asab', ncol = 2) + coord_fixed() + ggtitle('ATAC - Before Harmony')
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.atac", label = TRUE, split.by = 'hto.donor', ncol = 2) + coord_fixed() + ggtitle('ATAC - Before Harmony')
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.atac", label = TRUE, features = 'nCount_ATAC', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 5) - ATAC - Harmony integration ----
l2_teaseq_tcell_obj <- RunHarmony(
  object = l2_teaseq_tcell_obj,
  group.by.vars = 'hto.donor', # Integrating across donors
  assay.use = 'ATAC',
  reduction.use = 'lsi.atac',
  reduction.save = 'lsi.atac.harmony',
  project.dim = FALSE,
  dims.use = 2:30,
  seed = 26
)

# Inspect correlation between Harmony LSI components and ATAC read depth
FeatureScatter(l2_teaseq_tcell_obj, feature1 = 'lsiatacharmony_1', feature2 = 'nCount_ATAC')
FeatureScatter(l2_teaseq_tcell_obj, feature1 = 'lsiatacharmony_2', feature2 = 'nCount_ATAC')
ElbowPlot(l2_teaseq_tcell_obj, reduction = 'lsi.atac.harmony', ndims = 29)

# UMAP embedding and clustering of Harmonized ATAC reduction
l2_teaseq_tcell_obj <- FindNeighbors(l2_teaseq_tcell_obj, reduction = "lsi.atac.harmony", dims = 1:29, assay = 'ATAC')
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, dims = 1:29, reduction = "lsi.atac.harmony", reduction.name = "umap.atac.harmony", reduction.key = "atacharmonyUMAP_", assay = 'ATAC')
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, resolution = 0.5, cluster.name = "atac.tcell.subcluster.harmony", assay = 'ATAC')
Idents(l2_teaseq_tcell_obj) <- "atac.tcell.subcluster.harmony"

# Clustering visualization after Harmony
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.atac.harmony", label = TRUE) + ggtitle('ATAC - After Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.atac.harmony", label = TRUE, split.by = 'asab', ncol = 2) + coord_fixed() + ggtitle('ATAC - After Harmony')
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.atac.harmony", label = TRUE, split.by = 'hto.donor', ncol = 4) + coord_fixed() + ggtitle('ATAC - After Harmony')
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.atac.harmony", label = TRUE, features = 'nCount_ATAC', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.atac.harmony", label = TRUE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 6) RNA - dimensionality reduction before integration ----

# Standard Seurat workflow for DR and clustering
DefaultAssay(l2_teaseq_tcell_obj) <- "RNA"
l2_teaseq_tcell_obj <- NormalizeData(l2_teaseq_tcell_obj, assay = 'RNA')
l2_teaseq_tcell_obj <- FindVariableFeatures(l2_teaseq_tcell_obj, assay = 'RNA')
l2_teaseq_tcell_obj <- ScaleData(l2_teaseq_tcell_obj, assay = 'RNA')
l2_teaseq_tcell_obj <- RunPCA(l2_teaseq_tcell_obj, assay = 'RNA', reduction.name = 'pca.rna', reduction.key = 'rnaPC_')
ElbowPlot(l2_teaseq_tcell_obj, ndims = 50, reduction = 'pca.rna')
l2_teaseq_tcell_obj <- FindNeighbors(l2_teaseq_tcell_obj, assay = 'RNA', reduction = 'pca.rna', dims = 1:50)
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, assay = 'RNA', reduction = 'pca.rna', dims = 1:50, reduction.name = 'umap.rna')
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, assay = 'RNA', cluster.name = 'rna.tcell.subcluster', resolution = 1.7)
Idents(l2_teaseq_tcell_obj) <- 'rna.tcell.subcluster'

# Cluster visualization across biological and technical variables
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = TRUE, cols = tcell_cols_unspec) + ggtitle('RNA - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, group.by = 'asab', cols = tcell_cols_unspec) + ggtitle('RNA - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, split.by = 'asab', cols = tcell_cols_unspec) + ggtitle('RNA - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, group.by = 'hto.tissue', cols = tissue_cols) + ggtitle('RNA - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, split.by = 'hto.tissue', cols = tcell_cols_unspec) + ggtitle('RNA - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, group.by = 'CC.Phase', cols = tcell_cols_unspec) + ggtitle('RNA - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, split.by = 'CC.Phase', cols = tcell_cols_unspec) + ggtitle('RNA - Before Harmony') + coord_fixed()

# Feature exploration, including X and Y linked genes before integration
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, features = 'rna_GNG4', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, features = 'rna_CD4', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, features = 'rna_CD8A', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = FALSE, features = 'percent.ribo', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 7) RNA - Harmony integration ----

l2_teaseq_tcell_obj[["RNA"]] <- split(l2_teaseq_tcell_obj[["RNA"]], f = l2_teaseq_tcell_obj$hto.donor) # split RNA layers by donor
l2_teaseq_tcell_obj <- IntegrateLayers(object = l2_teaseq_tcell_obj, 
                                       method = HarmonyIntegration, 
                                       orig.reduction = "pca.rna", 
                                       new.reduction = "pca.rna.harmony",
                                       assay = 'RNA',
                                       verbose = TRUE,
                                       seed = 26)
l2_teaseq_tcell_obj <- FindNeighbors(l2_teaseq_tcell_obj, reduction = "pca.rna.harmony", dims = 1:50)
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, dims = 1:50, reduction = "pca.rna.harmony", reduction.name = "umap.rna.harmony", reduction.key = "rnaharmonyUMAP_")
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, resolution = 1, cluster.name = "rna.tcell.subcluster.harmony")
Idents(l2_teaseq_tcell_obj) <- "rna.tcell.subcluster.harmony"

# Visualization of integrated RNA embedding
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna.harmony", label = TRUE) + ggtitle('RNA - After Harmony Integration') + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_GNG4', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B'), split.by = 'hto.tissue') + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B'), split.by = 'hto.tissue') + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.rna.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B'), split.by = 'hto.tissue') + coord_fixed()

# Compare embeddings before versus after integration
Idents(l2_teaseq_tcell_obj) <- 'rna.tcell.subcluster'
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = TRUE) + ggtitle('RNA - Before Harmony Integration') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = TRUE, split.by = 'hto.tissue') + ggtitle('RNA - Before Harmony Integration') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = TRUE, split.by = 'asab') + ggtitle('RNA - Before Harmony Integration') + coord_fixed()
Idents(l2_teaseq_tcell_obj) <- "rna.tcell.subcluster.harmony"
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna.harmony", label = TRUE) + ggtitle('RNA - After Harmony Integration') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna.harmony", label = TRUE, split.by = 'hto.tissue') + ggtitle('RNA - After Harmony Integration') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.rna", label = TRUE, split.by = 'asab') + ggtitle('RNA - After Harmony Integration') + coord_fixed()

# 8) ADT - dimensionality reduction before integration ----

# Stardard Seurat workflow for ADT DR and clustering
DefaultAssay(l2_teaseq_tcell_obj) <- "ADT"
l2_teaseq_tcell_obj[["ADT"]] <- split(l2_teaseq_tcell_obj[["ADT"]], f = l2_teaseq_tcell_obj$hto.donor) # Split ADT layers by donor
l2_teaseq_tcell_obj <- NormalizeData(l2_teaseq_tcell_obj, assay= "ADT", normalization.method = "CLR", margin = 2)
l2_teaseq_tcell_obj <- ScaleData(l2_teaseq_tcell_obj, assay = "ADT")
l2_teaseq_tcell_obj <- FindVariableFeatures(l2_teaseq_tcell_obj, assay = "ADT")
l2_teaseq_tcell_obj <- RunPCA(l2_teaseq_tcell_obj, assay = "ADT", reduction.name = "pca.adt", reduction.key = "adtPC_")
ElbowPlot(l2_teaseq_tcell_obj, ndims = 50, reduction = "pca.adt")
l2_teaseq_tcell_obj <- FindNeighbors(l2_teaseq_tcell_obj, reduction = "pca.adt", dims = 1:50)
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, assay = "ADT", dims = 1:50, reduction = "pca.adt", reduction.name = "umap.adt", reduction.key = "adtUMAP_")
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, resolution = 0.5, cluster.name = "adt.tcell.subcluster", graph.name = "ADT_snn")
Idents(l2_teaseq_tcell_obj) <- "adt.tcell.subcluster"

# Visualization before integration
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.adt", label = TRUE) + coord_fixed() + ggtitle('ADT - Before Harmony')
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.adt", ncol = 2, label = TRUE, split.by = 'asab') + ggtitle('ADT - Before Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.adt", ncol = 4, label = TRUE, split.by = 'hto.donor') + ggtitle('ADT - Before Harmony') + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.adt", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.adt", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.adt", label = FALSE, features = 'rna_CD8A', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 9) ADT - harmony integration ----
l2_teaseq_tcell_obj <- IntegrateLayers(object = l2_teaseq_tcell_obj, 
                                method = HarmonyIntegration, 
                                orig.reduction = "pca.adt", 
                                new.reduction = "pca.adt.harmony",
                                assay = 'ADT',
                                verbose = TRUE,
                                seed = 26)
l2_teaseq_tcell_obj <- FindNeighbors(l2_teaseq_tcell_obj, reduction = "pca.adt.harmony", dims = 1:50, assay = 'ADT')
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, dims = 1:50, reduction = "pca.adt.harmony", reduction.name = "umap.adt.harmony", reduction.key = "adtharmonyUMAP_", assay = 'ADT')
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, resolution = 0.7, cluster.name = "adt.tcell.subcluster.harmony", assay = 'ADT')
Idents(l2_teaseq_tcell_obj) <- "adt.tcell.subcluster.harmony"

# Visualization after integration
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.adt.harmony", label = TRUE) + ggtitle('ADT - After Harmony') + coord_fixed()
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.adt.harmony", label = TRUE, split.by = 'asab', ncol = 2) + coord_fixed() + ggtitle('ADT - After Harmony')
DimPlot(l2_teaseq_tcell_obj, reduction = "umap.adt.harmony", label = TRUE, split.by = 'hto.donor', ncol = 2) + coord_fixed() + ggtitle('ADT - After Harmony')
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l2_teaseq_tcell_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_BCL6', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 10) 3WNN analysis using Harmony-corrected trimodal data ----

# Multimodal clustering and dimensionality reduction
l2_teaseq_tcell_obj <- FindMultiModalNeighbors(l2_teaseq_tcell_obj, 
                                        reduction.list = list("pca.rna.harmony", "lsi.atac.harmony", "pca.adt.harmony"), 
                                        dims.list = list(1:50, 1:29, 1:50))
l2_teaseq_tcell_obj <- RunUMAP(l2_teaseq_tcell_obj, nn.name = "weighted.nn", reduction.name = "umap.wnn.harmony", reduction.key = "harmonywnnUMAP_")
l2_teaseq_tcell_obj <- FindClusters(l2_teaseq_tcell_obj, graph.name = "wsnn", algorithm = 3, resolution = 1.08, verbose = TRUE, cluster.name = "wnn.tcell.subcluster.harmony")
Idents(l2_teaseq_tcell_obj) <- "wnn.tcell.subcluster.harmony"

# 11) Annotate TEAseq L2 3WNN clusters  ----

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

# Specify colors for each L2 T cell cluster
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

# 12) Add TF motif information to object ----

# Guided by Signac vignette - https://stuartlab.org/signac/articles/motif_vignette

# Filtering feature names - https://github.com/stuart-lab/signac/issues/780
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
gr <- granges(l2_teaseq_tcell_obj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
dim(l2_teaseq_tcell_obj)
l2_teaseq_tcell_obj[['ATAC']] <- subset(l2_teaseq_tcell_obj[["ATAC"]], features = feat.keep) # 227103  - 226955  = 148 peaks removed

# Get list of motif position frequency matrices from JASPAR database
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
pfm <- getMatrixSet(  # get TF motif PFM information
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information to object
l2_teaseq_tcell_obj <- AddMotifs(
  object = l2_teaseq_tcell_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# 13) Add ATAC ChromVAR assay to object ----

l2_teaseq_tcell_obj <- RunChromVAR(
  object = l2_teaseq_tcell_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = 'ATAC',
  new.assay.name = 'chromvar'
)

# 14) Add ATAC GeneActivity assay to L2 T cell object ----
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
tcell_geneactivity_mtx <- GeneActivity(l2_teaseq_tcell_obj,
                                       assay = 'ATAC',
                                       extend.upstream = 2000, # default
                                       extend.downstream = 0, # default
                                       features = rownames(l2_teaseq_tcell_obj@assays$RNA),
                                       biotypes = NULL, # include noncoding genes, otherwise default 'protein_coding'
                                       max.width = NULL, # do not filter based on gene length, otherwise default '5e+05'
                                       process_n = 2000, # decreasing runtime, default 2000
                                       gene.id = FALSE,
                                       verbose = TRUE
)
l2_teaseq_tcell_obj[['ACT']] <- CreateAssayObject(counts = tcell_geneactivity_mtx)
l2_teaseq_tcell_obj <- NormalizeData(
  object = l2_teaseq_tcell_obj,
  assay = 'ACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(l2_teaseq_tcell_obj$nCount_ACT)
)
l2_teaseq_tcell_obj <- ScaleData(l2_teaseq_tcell_obj, assay = 'ACT')

# 15) Add RNA SCENIC output to Seurat object - (refer to SCENIC folder for run information) ----
setwd('/filepath/l2_tcell_scenic_output_folder/')
scenicOptions <- readRDS("/filepath/l2_tcell_scenic_output_folder/scenicOptions4.rds") # after complete pipeline run, refer to SCENIC code folder
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC") # load AUC (regulon activity) matrix
AUCmatrix <- getAUC(regulonAUC)
l2_teaseq_tcell_obj[["SCENIC"]] <- CreateAssayObject(data = AUCmatrix)

# 16) Find TEAseq 3WNN L2 T Cell cluster markers across assays ----

# Join layers for differential expression testing
l2_teaseq_tcell_obj <- JoinLayers(l2_teaseq_tcell_obj, assay = 'RNA')
l2_teaseq_tcell_obj <- JoinLayers(l2_teaseq_tcell_obj, assay = 'ADT')
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'

# Find TEAseq 3WNN L2 T cell cluster markers for each assay using FindAllMarkers
# ATAC - peaks ('ATAC'), GeneActivity ('ACT', promoter region and gene body accessibility), and chromVAR (TF motif deviation Z-scores)
# RNA - transcripts ('RNA', log-normalized) and SCENIC (TF regulon AUC scores)
# ADT - antibody-derived tag epitope information (CLR-normalized)
# Default output column names were adjusted as detailed below to fit each assay and improve interpretability

# TEAseq 3WNN L2 T Cell Cluster Markers - ATAC (peaks)
DefaultAssay(l2_teaseq_tcell_obj) <- 'ATAC'
l2_tcell_wnn_clust_atac_peak_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ATAC') # using default parameters
dap_names <- l2_tcell_wnn_clust_atac_peak_markers$gene
dap_gr <- StringToGRanges(dap_names)
dap_closest_gene <- ClosestFeature(
  object = l2_teaseq_tcell_obj,
  regions    = dap_gr,
  annotation = Annotation(l2_teaseq_tcell_obj)
)
l2_tcell_wnn_clust_atac_peak_markers <- l2_tcell_wnn_clust_atac_peak_markers %>%
  mutate(
    closest_gene_symbol = dap_closest_gene$gene_name,
    closest_ensembl_id = dap_closest_gene$gene_id
  )
l2_tcell_wnn_clust_atac_peak_markers <- l2_tcell_wnn_clust_atac_peak_markers %>%
  rename(
    l2_tcell_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    atac_peak = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_tcell_wnn_clust_atac_peak_markers <- l2_tcell_wnn_clust_atac_peak_markers %>% arrange(l2_tcell_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tcell_wnn_clust_atac_peak_markers <- l2_tcell_wnn_clust_atac_peak_markers %>% relocate(c('l2_tcell_cluster','atac_peak','closest_gene_symbol','closest_ensembl_id','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tcell_wnn_clust_atac_peak_markers, '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_atac_peak_markers.rds')

# TEAseq 3WNN L2 T Cell Cluster Markers - ACT (ATAC-based Signac GeneActivity summary of gene body and promoter region accessibility)
DefaultAssay(l2_teaseq_tcell_obj) <- 'ACT'
l2_tcell_wnn_clust_atac_geneactivity_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ACT') # using default parameters
l2_tcell_wnn_clust_atac_geneactivity_markers <- l2_tcell_wnn_clust_atac_geneactivity_markers %>%
  rename(
    l2_tcell_cluster = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_tcell_wnn_clust_atac_geneactivity_markers <- l2_tcell_wnn_clust_atac_geneactivity_markers %>% arrange(l2_tcell_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tcell_wnn_clust_atac_geneactivity_markers <- l2_tcell_wnn_clust_atac_geneactivity_markers %>% relocate(c('l2_tcell_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tcell_wnn_clust_atac_geneactivity_markers, '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_atac_geneactivity_markers.rds')

# TEAseq 3WNN L2 T Cell Cluster Markers - chromVAR (ATAC-based TF motif deviation Z-scores)
l2_tcell_wnn_clust_atac_chromvar_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'chromvar', 
                                                           mean.fxn = rowMeans, # use row means function to find cluster differences as chromVAR Z-scores are not log-normalized like Seurat RNA data
                                                           fc.name = "mean_diff_z_score") # using default parameters
pfm <- getMatrixSet( 
  x = JASPAR2020, # get TF motif PFM information
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
id_to_tf <- sapply(pfm, function(x) name(x)) # function to get common TF name for each position weight matrix ID
names(id_to_tf) <- sapply(pfm, ID)
l2_tcell_wnn_clust_atac_chromvar_markers$tf_name <- id_to_tf[l2_tcell_wnn_clust_atac_chromvar_markers$gene]
l2_tcell_wnn_clust_atac_chromvar_markers <- l2_tcell_wnn_clust_atac_chromvar_markers %>%
  rename(
    l2_tcell_cluster = cluster,
    jaspar_pfm_id = gene,
    p_val_raw = p_val,
    pct_nonzero_deviation_in_clust = pct.1, # rename, as nonzero chromVAR Z-scores do not mean 'positive' as with RNA expression 
    pct_nonzero_deviation_in_others = pct.2
  )
l2_tcell_wnn_clust_atac_chromvar_markers <- l2_tcell_wnn_clust_atac_chromvar_markers %>% arrange(l2_tcell_cluster, desc(mean_diff_z_score), p_val_adj)
l2_tcell_wnn_clust_atac_chromvar_markers <- l2_tcell_wnn_clust_atac_chromvar_markers %>% relocate(c('l2_tcell_cluster','tf_name','jaspar_pfm_id','mean_diff_z_score','p_val_raw','p_val_adj'))
saveRDS(l2_tcell_wnn_clust_atac_chromvar_markers, '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_atac_chromvar_markers.rds')

# TEAseq 3WNN L2 T Cell Cluster Markers - RNA (standard Seurat log-normalized mRNA count data)
DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
l2_tcell_wnn_clust_rna_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'RNA') # using default parameters
l2_tcell_wnn_clust_rna_markers <- l2_tcell_wnn_clust_rna_markers %>%
  rename(
    l2_tcell_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_tcell_wnn_clust_rna_markers <- l2_tcell_wnn_clust_rna_markers %>% arrange(l2_tcell_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tcell_wnn_clust_rna_markers <- l2_tcell_wnn_clust_rna_markers %>% relocate(c('l2_tcell_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tcell_wnn_clust_rna_markers, '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_rna_markers.rds')

# TEAseq 3WNN L2 T Cell Cluster Markers - SCENIC (RNA-based TF regulon inference)
DefaultAssay(l2_teaseq_tcell_obj) <- "SCENIC"
l2_tcell_wnn_clust_rna_scenic_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'SCENIC', 
                                                        mean.fxn = rowMeans, # SCENIC regulon AUC scores are not log-transformed
                                                        fc.name = "mean_diff_AUC",
                                                        min.pct = 0, logfc.threshold = 0) # adjust default filters as the mean difference in AUC may be lower than 0.1 
l2_tcell_wnn_clust_rna_scenic_markers <- l2_tcell_wnn_clust_rna_scenic_markers %>%
  rename(
    l2_tcell_cluster = cluster,
    scenic_regulon = gene,
    p_val_raw = p_val,
    pct_nonzero_auc_in_clust = pct.1,
    pct_nonzero_auc_others = pct.2
  )
l2_tcell_wnn_clust_rna_scenic_markers <- l2_tcell_wnn_clust_rna_scenic_markers %>% arrange(l2_tcell_cluster, desc(mean_diff_AUC), p_val_adj)
l2_tcell_wnn_clust_rna_scenic_markers <- l2_tcell_wnn_clust_rna_scenic_markers %>% relocate(c('l2_tcell_cluster','scenic_regulon','mean_diff_AUC','p_val_raw','p_val_adj'))
saveRDS(l2_tcell_wnn_clust_rna_scenic_markers, '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_rna_scenic_markers.rds')

# TEAseq 3WNN L2 T Cell Cluster Markers - ADT (CLR-normalized count data)
DefaultAssay(l2_teaseq_tcell_obj) <- 'ADT'
l2_tcell_wnn_clust_adt_markers <- FindAllMarkers(l2_teaseq_tcell_obj, assay = 'ADT') # using default parameters
l2_tcell_wnn_clust_adt_markers <- l2_tcell_wnn_clust_adt_markers %>%
  rename(
    l2_tcell_cluster = cluster,
    adt_epitope = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l2_tcell_wnn_clust_adt_markers <- l2_tcell_wnn_clust_adt_markers %>% arrange(l2_tcell_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l2_tcell_wnn_clust_adt_markers$adt_epitope <- gsub("^adt-", "", l2_tcell_wnn_clust_adt_markers$adt_epitope)
l2_tcell_wnn_clust_adt_markers <- l2_tcell_wnn_clust_adt_markers %>% relocate(c('l2_tcell_cluster','adt_epitope','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l2_tcell_wnn_clust_adt_markers, '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_adt_markers.rds')

# Batch import .rds marker files and export as single multi-sheet XLSX file
clust_marker_df_dir   <- '/filepath/step14_l2_teaseq_tcell_analysis/'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = '/filepath/step14_l2_teaseq_tcell_analysis/l2_tcell_wnn_clust_markers_all_assays.xlsx')  


# 17) Save L2 object and proceed to L3 subclustering of Tfh-like cells ----
saveRDS(l2_teaseq_tcell_obj, '/filepath/step14_l2_teaseq_tcell_analysis/l2_teaseq_tcell_obj.rds')