# Step 13 - Level 1 TEAseq Analysis - All mononuclear cells in tonsil and peripheral blood samples
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
setwd("/filepath/step13_l1_teaseq_analysis/")

# Colors assigned to Level 1 3WNN clusters
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

# Unspecified colors
l1_cols_unspec <- c(
  "#8ab1cc", "#657ab0", "#94a890", "#add5df", "#e49a78",
  "#d4b5e4", "#c49980", "#728782", "#c9a3c1", "#d66d97",
  "#b6b7ba", "#e0c1b4", "#e1d4c7", "#d699ac", "#dfcc78",
  "#bad4f4", "#eaebeb", "#f4a6c7", "#e17794", "#b3f5dd",
  "#535c81", "#50664f", "#4a3c5d", "#72d1b4", "#f7f1e0",
  "#c9744d", "#b37cbf"
)

# Numbers assigned colors 0-26
l1_cols_numb <- c(
  "0"  = "#8ab1cc", "1"  = "#657ab0", "2"  = "#94a890",
  "3"  = "#add5df", "4"  = "#e49a78", "5"  = "#d4b5e4",
  "6"  = "#c49980", "7"  = "#728782", "8"  = "#c9a3c1",
  "9"  = "#d66d97", "10" = "#b6b7ba", "11" = "#e0c1b4",
  "12" = "#e1d4c7", "13" = "#d699ac", "14" = "#dfcc78",
  "15" = "#bad4f4", "16" = "#eaebeb", "17" = "#f4a6c7",
  "18" = "#e17794", "19" = "#b3f5dd", "20" = "#535c81",
  "21" = "#50664f", "22" = "#4a3c5d", "23" = "#72d1b4",
  "24" = "#f7f1e0", "25" = "#c9744d", "26" = "#b37cbf"
)

# 2) Import Seurat object from TEAseq Data Preprocessing Step 12 (all tonsil and peripheral blood mononuclear cells - post multiplet filtering, QC filtering, peak matrix unification, GEM well merging, and cell-cycle scoring) ----
l1_teaseq_obj <- readRDS("/filepath/step12_rna_cc_scoring/merged.cc.scored.rds")

# 3) ATAC - dimensionality reduction before integration ----
DefaultAssay(l1_teaseq_obj) <- "ATAC"
l1_teaseq_obj <- RunTFIDF(l1_teaseq_obj, assay = 'ATAC')
l1_teaseq_obj <- FindTopFeatures(l1_teaseq_obj, min.cutoff = 'q0', assay = 'ATAC')
l1_teaseq_obj <- RunSVD(l1_teaseq_obj, reduction.name = "lsi.atac", reduction.key = "atacLSI_", assay = 'ATAC')

# Selection of components for integration
ElbowPlot(l1_teaseq_obj, ndims = 50, reduction = "lsi.atac")
FeatureScatter(l1_teaseq_obj, feature1 = 'atacLSI_1', feature2 = 'nCount_ATAC') # First ATAC LSI component is strongly correlated with read depth as expected
FeatureScatter(l1_teaseq_obj, feature1 = 'atacLSI_2', feature2 = 'nCount_ATAC')
FeatureScatter(l1_teaseq_obj, feature1 = 'atacLSI_3', feature2 = 'nCount_ATAC')
FeatureScatter(l1_teaseq_obj, feature1 = 'atacLSI_4', feature2 = 'nCount_ATAC')

# Clustering and UMAP embedding
l1_teaseq_obj <- FindNeighbors(l1_teaseq_obj, reduction= "lsi.atac", dims = 2:50, assay = 'ATAC')
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, resolution = 1, cluster.name = "atac.bulk.cluster", graph.name = "ATAC_snn")
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, reduction = 'lsi.atac', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_", assay = 'ATAC')

# Visualizing technical and biological variables before integration
Idents(l1_teaseq_obj) <- 'atac.bulk.cluster'
DimPlot(l1_teaseq_obj, reduction = "umap.atac", label = TRUE) + ggtitle('ATAC - Before Harmony') + coord_fixed()
DimPlot(l1_teaseq_obj, reduction = "umap.atac", label = TRUE, split.by = 'asab', ncol = 2) + coord_fixed() + ggtitle('ATAC - Before Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.atac", label = TRUE, split.by = 'hto.donor', ncol = 2) + coord_fixed() + ggtitle('ATAC - Before Harmony')
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac", label = TRUE, features = 'nCount_ATAC', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac", label = TRUE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac", label = FALSE, features = 'percent.ribo', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac", label = FALSE, features = 'percent.mt', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 4) ATAC - Harmony integration ----
DefaultAssay(l1_teaseq_obj) <- 'ATAC'
l1_teaseq_obj <- RunHarmony(
  object = l1_teaseq_obj,
  group.by.vars = 'hto.donor', # Integrating across donors
  assay.use = 'ATAC',
  reduction.use = 'lsi.atac',
  reduction.save = 'lsi.atac.harmony',
  project.dim = FALSE,
  dims.use = 2:50,
  seed = 26
)

# Assessing Harmonized LSI components after integration
FeatureScatter(l1_teaseq_obj, feature1 = 'lsiatacharmony_1', feature2 = 'nCount_ATAC')
FeatureScatter(l1_teaseq_obj, feature1 = 'lsiatacharmony_2', feature2 = 'nCount_ATAC')
ElbowPlot(l1_teaseq_obj, reduction = 'lsi.atac.harmony', ndims = 49)

# Clustering and UMAP embedding after integration
l1_teaseq_obj <- FindNeighbors(l1_teaseq_obj, reduction = "lsi.atac.harmony", dims = 1:49, assay = 'ATAC', graph.name = c('ATAC_harmony_nn','ATAC_harmony_snn'))
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, dims = 1:49, reduction = "lsi.atac.harmony", reduction.name = "umap.atac.harmony", reduction.key = "atacharmonyUMAP_", assay = 'ATAC')
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, resolution = 0.45, cluster.name = "atac.bulk.cluster.harmony", assay = 'ATAC', graph.name = 'ATAC_harmony_snn')
Idents(l1_teaseq_obj) <- "atac.bulk.cluster.harmony"
DimPlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = TRUE) + ggtitle('ATAC - After Harmony') + coord_fixed()

# 5) ATAC - Cluster Annotation ----

# Level 1 ATAC Cluster Annotation (Coarse - 15 Lineages)

bulk_atac_cluster_names <- c(
  "0" = "CD4 Tn 1",
  "1" = "CD4 Tn 2",
  "2" = "Tfh-like",
  "3" = "NBC 1",
  "4" = "CD4 Tm",
  "5" = "CD8 Tm/Inn",
  "6" = "CD4 Tn 3",
  "7" = "MBC",
  "8" = "CD8 Tn",
  "9" = "NK",
  "10" = "GCB",
  "11" = "Myl",
  "12" = "NBC 2",
  "13" = "Treg",
  "14" = "ASC"
)

atac_cluster_colors <- c(
  "CD4 Tn 1" = "#8ab1cc",       
  "NBC 1" = "#90a6bf",        
  "CD4 Tm" = "#94a890",      
  "CD8 Tm/Inn" = "#add5df",
  "NBC 2" = "#b3f5dd",
  "Tfh-like" = "#e49a78",         
  "MBC" = "#6e485a",         
  "CD8 Tn" = "#50664f",      
  "NK" = "#723d56",          
  "Treg" = "#c9a3c1",        
  "GCB" = "#d66d97",         
  "Myl" = "#535c81",         
  "Ribo" = "#e0c1b4",        
  "CD4 Tn 2" = "#e1d4c7",    
  "CD4 Tn 3" = "#7b7f9c",  
  "ASC" = "#72d1b4"          
)

# Rename Level 1 coarse 15 ATAC cluster numbers to assigned names
l1_teaseq_obj$bulk.atac.annot <- l1_teaseq_obj$atac.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'bulk.atac.annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, bulk_atac_cluster_names)
l1_teaseq_obj$bulk.atac.annot <- Idents(l1_teaseq_obj)

# Verify renamed clusters
table(Idents(l1_teaseq_obj))

# Visualizing biological and technical variables after Harmony integration
DimPlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = TRUE, split.by = 'asab', ncol = 2) + coord_fixed() + ggtitle('ATAC - After Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = TRUE, split.by = 'hto.donor', ncol = 4) + coord_fixed() + ggtitle('ATAC - After Harmony')
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = TRUE, features = 'nCount_ATAC', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = TRUE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'percent.ribo', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.atac.harmony", label = FALSE, features = 'percent.mt', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 6) RNA - dimensionality reduction before integration ----

DefaultAssay(l1_teaseq_obj) <- "RNA"
l1_teaseq_obj <- NormalizeData(l1_teaseq_obj, assay = 'RNA')
l1_teaseq_obj <- FindVariableFeatures(l1_teaseq_obj, assay = 'RNA')
l1_teaseq_obj <- ScaleData(l1_teaseq_obj, assay = 'RNA')
l1_teaseq_obj <- RunPCA(l1_teaseq_obj, assay = 'RNA', reduction.name = 'pca.rna', reduction.key = 'rnaPC_')
ElbowPlot(l1_teaseq_obj, ndims = 50, reduction = 'pca.rna')
l1_teaseq_obj <- FindNeighbors(l1_teaseq_obj, assay = 'RNA', reduction = 'pca.rna', dims = 1:50)
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, assay = 'RNA', reduction = 'pca.rna', dims = 1:50, reduction.name = 'umap.rna')
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, assay = 'RNA', cluster.name = 'rna.bulk.cluster', resolution = 0.45, graph.name = 'RNA_snn')
Idents(l1_teaseq_obj) <- 'rna.bulk.cluster'
DimPlot(l1_teaseq_obj, reduction = "umap.rna", label = TRUE) + ggtitle('RNA - Before Harmony') + coord_fixed()

# 7) RNA - Harmony integration ----

DefaultAssay(l1_teaseq_obj) <- "RNA"
DimPlot(l1_teaseq_obj, reduction = 'umap.rna', label = TRUE) + coord_fixed()
l1_teaseq_obj[["RNA"]] <- split(l1_teaseq_obj[["RNA"]], f = l1_teaseq_obj$hto.donor) # split RNA layers by donor
l1_teaseq_obj <- IntegrateLayers(object = l1_teaseq_obj, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca.rna", 
                                 new.reduction = "pca.rna.harmony",
                                 assay = 'RNA',
                                 verbose = TRUE,
                                 seed = 26)
l1_teaseq_obj <- FindNeighbors(l1_teaseq_obj, assay = 'RNA', reduction = "pca.rna.harmony", dims = 1:50, graph.name = c('RNA_harmony_nn','RNA_harmony_snn'))
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, assay = 'RNA', dims = 1:50, reduction = "pca.rna.harmony", reduction.name = "umap.rna.harmony", reduction.key = "rnaharmonyUMAP_")
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, assay = 'RNA', resolution = 0.45, cluster.name = "rna.bulk.cluster.harmony", graph.name = 'RNA_harmony_snn')
Idents(l1_teaseq_obj) <- "rna.bulk.cluster.harmony"
DimPlot(l1_teaseq_obj, reduction = "umap.rna.harmony", label = TRUE) + ggtitle('RNA - After Harmony Integration') + coord_fixed()

# Compare embeddings
Idents(l1_teaseq_obj) <- 'rna.bulk.cluster'
DimPlot(l1_teaseq_obj, reduction = "umap.rna", label = TRUE) + ggtitle('RNA - Before Harmony Integration') + coord_fixed()
DimPlot(l1_teaseq_obj, reduction = "umap.rna", label = TRUE, split.by = 'hto.tissue') + ggtitle('RNA - Before Harmony Integration') + coord_fixed()
DimPlot(l1_teaseq_obj, reduction = "umap.rna", label = TRUE, split.by = 'asab') + ggtitle('RNA - Before Harmony Integration') + coord_fixed()
Idents(l1_teaseq_obj) <- "rna.bulk.cluster.harmony"
DimPlot(l1_teaseq_obj, reduction = "umap.rna.harmony", label = TRUE) + ggtitle('RNA - After Harmony Integration') + coord_fixed()

# 8) RNA - Cluster Annotation ----

l1_teaseq_obj <- JoinLayers(l1_teaseq_obj, assay = 'RNA')
Idents(l1_teaseq_obj) <- 'rna.bulk.cluster.harmony'
table(l1_teaseq_obj$rna.bulk.cluster.harmony)

# Level 1 RNA Cluster Annotation (Coarse - 15 Lineages)

bulk_rna_cluster_names <- c(
  "0" = "CD4 Tn",
  "1" = "CD4 Tm",
  "2" = "NBC",
  "3" = "CD8 Tm/Inn 1",
  "4" = "Tfh-like",
  "5" = "MBC",
  "6" = "CD8 Tn",
  "7" = "Treg",
  "8" = "GCB",
  "9" = "NK",
  "10" = "CD8 Tm/Inn 2",
  "11" = "Myl",
  "12" = "Myl CD16",
  "13" = "Ribo",
  "14" = "ASC"
)

rna_cluster_colors <- c(
  "CD4 Tn" = "#8ab1cc",       
  "NBC" = "#90a6bf",        
  "CD4 Tm" = "#94a890",      
  "CD8 Tm/Inn 1" = "#add5df",
  "CD8 Tm/Inn 2" = "#b3f5dd",
  "Tfh-like" = "#e49a78",         
  "MBC" = "#6e485a",         
  "CD8 Tn" = "#50664f",      
  "NK" = "#723d56",          
  "Treg" = "#c9a3c1",        
  "GCB" = "#d66d97",         
  "Myl" = "#535c81",         
  "Ribo" = "#e0c1b4",        
  "Myl CD16" = "#e1d4c7",    
  "CD8 Tunc" = "#7b7f9c",  
  "ASC" = "#72d1b4"          
)

# Rename Level 1 coarse 15 RNA cluster numbers to assigned names
l1_teaseq_obj$bulk.rna.annot <- l1_teaseq_obj$rna.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'bulk.rna.annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, bulk_rna_cluster_names)
l1_teaseq_obj$bulk.rna.annot <- Idents(l1_teaseq_obj)
table(Idents(l1_teaseq_obj))

# 9) ADT - dimensionality reduction before integration ----

DefaultAssay(l1_teaseq_obj) <- "ADT"
l1_teaseq_obj[["ADT"]] <- split(l1_teaseq_obj[["ADT"]], f = l1_teaseq_obj$hto.donor) # Split ADT layers by donor
l1_teaseq_obj <- NormalizeData(l1_teaseq_obj, assay= "ADT", normalization.method = "CLR", margin = 2)
l1_teaseq_obj <- ScaleData(l1_teaseq_obj, assay = "ADT")
l1_teaseq_obj <- FindVariableFeatures(l1_teaseq_obj, assay = "ADT")
l1_teaseq_obj <- RunPCA(l1_teaseq_obj, assay = "ADT", reduction.name = "pca.adt", reduction.key = "adtPC_")
ElbowPlot(l1_teaseq_obj, ndims = 50, reduction = "pca.adt")
l1_teaseq_obj <- FindNeighbors(l1_teaseq_obj, reduction = "pca.adt", dims = 1:50)
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, assay = "ADT", dims = 1:50, reduction = "pca.adt", reduction.name = "umap.adt", reduction.key = "adtUMAP_")
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, resolution = 0.84, cluster.name = "adt.bulk.cluster", graph.name = "ADT_snn", assay = 'ADT')
Idents(l1_teaseq_obj) <- "adt.bulk.cluster"
DimPlot(l1_teaseq_obj, reduction = "umap.adt", label = TRUE) + coord_fixed() + ggtitle('ADT - Before Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.adt", ncol = 2, label = TRUE, split.by = 'asab') + ggtitle('ADT - Before Harmony') + coord_fixed()
DimPlot(l1_teaseq_obj, reduction = "umap.adt", ncol = 4, label = TRUE, split.by = 'hto.donor') + ggtitle('ADT - Before Harmony') + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.adt", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.adt", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 10) ADT - Harmony integration ----
l1_teaseq_obj <- IntegrateLayers(object = l1_teaseq_obj, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca.adt", 
                                 new.reduction = "pca.adt.harmony",
                                 assay = 'ADT',
                                 verbose = TRUE,
                                 seed = 26)
l1_teaseq_obj <- FindNeighbors(l1_teaseq_obj, reduction = "pca.adt.harmony", dims = 1:50, assay = 'ADT',  graph.name = c('ADT_harmony_nn','ADT_harmony_snn'))
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, dims = 1:50, reduction = "pca.adt.harmony", reduction.name = "umap.adt.harmony", reduction.key = "adtharmonyUMAP_", assay = 'ADT')
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, resolution = 0.25, cluster.name = "adt.bulk.cluster.harmony", assay = 'ADT', graph.name = 'ADT_harmony_snn')
Idents(l1_teaseq_obj) <- "adt.bulk.cluster.harmony"
DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = TRUE) + ggtitle('ADT - After Harmony') + coord_fixed()

# 11) ADT - Cluster Annotation ----

l1_teaseq_obj <- JoinLayers(l1_teaseq_obj, assay = 'ADT')
Idents(l1_teaseq_obj) <- 'adt.bulk.cluster.harmony'

# Level 1 ADT Cluster Annotation (Coarse - 15 Lineages)

bulk_adt_cluster_names <- c(
  "0" = "CD4 Tn 1",
  "1" = "NBC",
  "2" = "CD4 Tm",
  "3" = "Tfh-like",
  "4" = "CD8 Tm/Inn 1",
  "5" = "MBC/GCB",
  "6" = "CD8 Tn",
  "7" = "NK",
  "8" = "CD8 Tm/Inn 2",
  "9" = "CD4 Tn 2",
  "10" = "DC",
  "11" = "Mono CD16",
  "12" = "CD8 Tm/Inn 3",
  "13" = "Mono CD14",
  "14" = "MAIT"
)

adt_cluster_colors <- c(
  "CD4 Tn 1" = "#8ab1cc",       
  "NBC" = "#90a6bf",        
  "CD4 Tm" = "#94a890",      
  "CD8 Tm/Inn 1" = "#add5df",
  "CD8 Tm/Inn 2" = "#b3f5dd",
  "Tfh-like" = "#e49a78",         
  "MBC/GCB" = "#6e485a",         
  "CD8 Tn" = "#50664f",      
  "NK" = "#723d56",          
  "CD8 Tm/Inn 3" = "#c9a3c1",        
  "Myl 1" = "#d66d97",         
  "Myl 2" = "#535c81",         
  "MAIT N" = "#e0c1b4",        
  "Myl CD16" = "#e1d4c7",    
  "CD4 Tn 2" = "#7b7f9c"
)

# Rename Level 1 coarse 15 ADT cluster numbers to assigned names
l1_teaseq_obj$bulk.adt.annot <- l1_teaseq_obj$adt.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'bulk.adt.annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, bulk_adt_cluster_names)
l1_teaseq_obj$bulk.adt.annot <- Idents(l1_teaseq_obj)
table(Idents(l1_teaseq_obj))

# Visualize clusters across biological and technical variables
DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = TRUE, split.by = 'asab', ncol = 2) + coord_fixed() + ggtitle('ADT - After Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = TRUE, split.by = 'hto.donor', ncol = 4) + coord_fixed() + ggtitle('ADT - After Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = FALSE, split.by = 'hto.sort', ncol = 4) + coord_fixed() + ggtitle('ADT - After Harmony')
FeaturePlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.adt.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 12) 3WNN analysis using Harmony corrected modalities ----

# Trimodal clustering and UMAP embedding
l1_teaseq_obj <- FindMultiModalNeighbors(l1_teaseq_obj, 
                                         reduction.list = list("pca.rna.harmony", "lsi.atac.harmony", "pca.adt.harmony"), 
                                         dims.list = list(1:50, 1:49, 1:50))
l1_teaseq_obj <- RunUMAP(l1_teaseq_obj, nn.name = "weighted.nn", reduction.name = "umap.wnn.harmony", reduction.key = "harmonywnnUMAP_")
l1_teaseq_obj <- FindClusters(l1_teaseq_obj, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = TRUE, cluster.name = "wnn.bulk.cluster.harmony")
Idents(l1_teaseq_obj) <- "wnn.bulk.cluster.harmony"

# Visualize clusters from 3WNN analysis of Harmony integrated trimodal data
Idents(l1_teaseq_obj) <- "l1_wnn_annot"
DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = TRUE, pt.size = 0.5) + coord_fixed() + ggtitle('3WNN - After Harmony')

# Visualize clusters across biological and technical variables
DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = TRUE, pt.size = 0.5, split.by = 'hto.tissue', ncol = 2) + coord_fixed() + ggtitle('3WNN - After Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = TRUE, pt.size = 0.5, split.by = 'asab') + coord_fixed() + ggtitle('3WNN - After Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = TRUE, pt.size = 0.5, split.by = 'hto.donor', ncol = 4) + coord_fixed() + ggtitle('3WNN - After Harmony')
DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = FALSE, group.by = 'hto.tissue', cols = c('steelblue2','orange2'), pt.size = 0.75, alpha = 0.75, order = FALSE) + coord_fixed() + ggtitle('') + NoLegend()
FeaturePlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'rna_UTY', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()
FeaturePlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = FALSE, features = 'rna_XIST', pt.size = 1, order = TRUE, cols = c('lightgrey','#7F1E2B')) + coord_fixed()

# 13) Annotate 3WNN L1 clusters ----

# 3WNN L1 Annotation
Idents(l1_teaseq_obj) <- "wnn.bulk.cluster.harmony"
l1_wnn_clust_names <- c(
  "0" = "CD4 Tn",
  "1" = "NBC",
  "2" = "CD4 Tm",
  "3" = "CD8 Tm/Inn",
  "4" = "Tfh-like",
  "5" = "MBC",
  "6" = "CD8 Tn",
  "7" = "NK",
  "8" = "Treg",
  "9" = "GCB",
  "10" = "Myl",
  "11" = "CD4 Tn Ribo",
  "12" = "Myl CD16",
  "13" = "CD8 UTC",
  "14" = "ASC"
)

# Rename 15 coarse 3WNN cluster numbers to assigned names
l1_teaseq_obj$l1_wnn_annot <- l1_teaseq_obj$wnn.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_wnn_clust_names)
l1_teaseq_obj$l1_wnn_annot <- Idents(l1_teaseq_obj)

# 14) Add TF motif information to object ----

# Guided by Signac vignette - https://stuartlab.org/signac/articles/motif_vignette

# Filtering feature names - https://github.com/stuart-lab/signac/issues/780
DefaultAssay(l1_teaseq_obj) <- 'ATAC'
gr <- granges(l1_teaseq_obj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
l1_teaseq_obj[['ATAC']] <- subset(l1_teaseq_obj[["ATAC"]], features = feat.keep) # 148 Features filtered - 227103-226955

# Get list of motif position frequency matrices from JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information to object
l1_teaseq_obj <- AddMotifs(
  object = l1_teaseq_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# 15) Add ATAC ChromVAR assay to object ----

register(MulticoreParam(28)) 
l1_teaseq_obj <- RunChromVAR(
  object = l1_teaseq_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = 'ATAC',
  new.assay.name = 'chromvar'
)

# 16) Add ATAC GeneActivity assay to object ----

DefaultAssay(l1_teaseq_obj) <- 'ATAC'
l1_teaseq_geneactivity_mtx <- GeneActivity(l1_teaseq_obj,
                                           assay = 'ATAC',
                                           extend.upstream = 2000, # default
                                           extend.downstream = 0, # default
                                           biotypes = NULL, # include noncoding genes, otherwise default 'protein_coding'
                                           max.width = NULL, # do not filter based on gene length, otherwise default '5e+05'
                                           process_n = 2000, # decreasing runtime, default 2000
                                           gene.id = FALSE,
                                           verbose = TRUE,
                                           features = rownames(l1_teaseq_obj@assays$RNA)
)
l1_teaseq_obj[['ACT']] <- CreateAssayObject(counts = l1_teaseq_geneactivity_mtx)
l1_teaseq_obj <- NormalizeData(
  object = l1_teaseq_obj,
  assay = 'ACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(l1_teaseq_obj$nCount_ACT) # following Signac vignette
)
l1_teaseq_obj <- ScaleData(l1_teaseq_obj, assay = 'ACT')

# 17) Add RNA SCENIC output to Seurat object - (refer to SCENIC folder for run information) ----
setwd('/filepath/l1_scenic_output_folder/')
scenicOptions <- readRDS("/filepath/l1_scenic_output_folder/scenicOptions4.rds") # after complete pipeline run, refer to SCENIC code folder
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC") # load AUC (regulon activity) matrix
AUCmatrix <- getAUC(regulonAUC)
l1_teaseq_obj[["SCENIC"]] <- CreateAssayObject(data = AUCmatrix)
DefaultAssay(l1_teaseq_obj) <- "SCENIC"

# 18) Find L1 TEAseq 3WNN cluster markers across assays ----

# Set cluster identity to 3WNN L1 cluster name
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'

# Join RNA and ADT layers to run FindAllMarkers - other assays do not have layers split by donor
l1_teaseq_obj <- JoinLayers(l1_teaseq_obj, assay = 'RNA')
l1_teaseq_obj <- JoinLayers(l1_teaseq_obj, assay = 'ADT')

# Find TEAseq 3WNN L1 global cluster markers for each assay using FindAllMarkers, adjusted as detailed below to fit the data type of each assay
# ATAC - peaks ('ATAC'), GeneActivity ('ACT', promoter region and gene body accessibility), and chromVAR (TF motif deviation scores)
# RNA - transcripts ('RNA', log-normalized) and SCENIC (TF regulon AUC scores)
# ADT - antibody-derived tag epitope information (CLR-normalized)

# TEAseq 3WNN L1 Cluster Markers - ATAC (peaks)
DefaultAssay(l1_teaseq_obj) <- 'ATAC'
l1_wnn_clust_atac_peak_markers <- FindAllMarkers(l1_teaseq_obj, assay = 'ATAC') # using default parameters
dap_names <- l1_wnn_clust_atac_peak_markers$gene
dap_gr <- StringToGRanges(dap_names)
Annotation(l1_teaseq_obj) <- gene_annot_hg38
dap_closest_gene <- ClosestFeature(
  object = l1_teaseq_obj,
  regions    = dap_gr,
  annotation = Annotation(l1_teaseq_obj)
)
l1_wnn_clust_atac_peak_markers <- l1_wnn_clust_atac_peak_markers %>%
  mutate(
    closest_gene_symbol = dap_closest_gene$gene_name,
    closest_ensembl_id = dap_closest_gene$gene_id
  )
l1_wnn_clust_atac_peak_markers <- l1_wnn_clust_atac_peak_markers %>%
  rename(
    l1_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    atac_peak = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_wnn_clust_atac_peak_markers <- l1_wnn_clust_atac_peak_markers %>% arrange(l1_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_wnn_clust_atac_peak_markers <- l1_wnn_clust_atac_peak_markers %>% relocate(c('l1_cluster','atac_peak','closest_gene_symbol','closest_ensembl_id','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_wnn_clust_atac_peak_markers, '/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_atac_peak_markers.rds')

# TEAseq 3WNN L1 Cluster Markers - ACT (ATAC-based Signac GeneActivity summary of gene body and promoter region accessibility)
DefaultAssay(l1_teaseq_obj) <- 'ACT'
l1_wnn_clust_atac_geneactivity_markers <- FindAllMarkers(l1_teaseq_obj, assay = 'ACT') # using default parameters
l1_wnn_clust_atac_geneactivity_markers <- l1_wnn_clust_atac_geneactivity_markers %>%
  rename(
    l1_cluster = cluster,
    gene_symbol = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_wnn_clust_atac_geneactivity_markers <- l1_wnn_clust_atac_geneactivity_markers %>% arrange(l1_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_wnn_clust_atac_geneactivity_markers <- l1_wnn_clust_atac_geneactivity_markers %>% relocate(c('l1_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_wnn_clust_atac_geneactivity_markers, '/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_atac_geneactivity_markers.rds')

# TEAseq 3WNN L1 Cluster Markers - chromVAR (ATAC-based TF motif deviation Z-scores)
l1_wnn_clust_atac_chromvar_markers <- FindAllMarkers(l1_teaseq_obj, assay = 'chromvar', 
                                                     mean.fxn = rowMeans, # use row means function to find cluster differences as chromVAR Z-scores are not log-normalized like Seurat RNA data
                                                     fc.name = "mean_diff_z_score") # using default parameters
pfm <- getMatrixSet( # get TF motif PFM information
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
id_to_tf <- sapply(pfm, function(x) name(x)) # function to get common TF name for each position weight matrix ID
names(id_to_tf) <- sapply(pfm, ID)
l1_wnn_clust_atac_chromvar_markers$tf_name <- id_to_tf[l1_wnn_clust_atac_chromvar_markers$gene]
l1_wnn_clust_atac_chromvar_markers <- l1_wnn_clust_atac_chromvar_markers %>%
  rename(
    l1_cluster = cluster,
    jaspar_pfm_id = gene,
    p_val_raw = p_val,
    pct_nonzero_deviation_in_clust = pct.1, # rename, as nonzero chromVAR Z-scores do not mean 'positive' as with RNA expression 
    pct_nonzero_deviation_in_others = pct.2
  )
l1_wnn_clust_atac_chromvar_markers <- l1_wnn_clust_atac_chromvar_markers %>% arrange(l1_cluster, desc(mean_diff_z_score), p_val_adj)
l1_wnn_clust_atac_chromvar_markers <- l1_wnn_clust_atac_chromvar_markers %>% relocate(c('l1_cluster','tf_name','jaspar_pfm_id','mean_diff_z_score','p_val_raw','p_val_adj'))
saveRDS(l1_wnn_clust_atac_chromvar_markers, '/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_atac_chromvar_markers.rds')

# TEAseq 3WNN L1 Cluster Markers - RNA (standard Seurat log-normalized mRNA count data)
DefaultAssay(l1_teaseq_obj) <- 'RNA'
l1_wnn_clust_rna_markers <- FindAllMarkers(l1_teaseq_obj, assay = 'RNA') # using default parameters
<- l1_wnn_clust_rna_markers %>%
  rename(
    l1_cluster = cluster,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_wnn_clust_rna_markers <- l1_wnn_clust_rna_markers %>% arrange(l1_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_wnn_clust_rna_markers <- l1_wnn_clust_rna_markers %>% relocate(c('l1_cluster','gene_symbol','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_wnn_clust_rna_markers, '/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_rna_markers.rds')

# SCENIC Markers (code below computes differential regulon activity between clusters - refer to SCENIC folder for run information)
l1_wnn_clust_rna_scenic_markers <- FindAllMarkers(l1_teaseq_obj, assay = 'SCENIC', 
                                                  mean.fxn = rowMeans, # SCENIC regulon AUC scores are not log-transformed
                                                  fc.name = "mean_diff_AUC",
                                                  min.pct = 0, logfc.threshold = 0) # adjust default filters as the mean difference in AUC may be lower than 0.1 
l1_wnn_clust_rna_scenic_markers <- l1_wnn_clust_rna_scenic_markers %>%
  rename(
    l1_cluster = cluster,
    scenic_regulon = gene,
    p_val_raw = p_val,
    pct_nonzero_auc_in_clust = pct.1,
    pct_nonzero_auc_others = pct.2
  )
l1_wnn_clust_rna_scenic_markers <- l1_wnn_clust_rna_scenic_markers %>% arrange(l1_cluster, desc(mean_diff_AUC), p_val_adj)
l1_wnn_clust_rna_scenic_markers <- l1_wnn_clust_rna_scenic_markers %>% relocate(c('l1_cluster','scenic_regulon','mean_diff_AUC','p_val_raw','p_val_adj'))
saveRDS(l1_wnn_clust_rna_scenic_markers, '/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_rna_scenic_markers.rds')

# TEAseq 3WNN L1 Cluster Markers - ADT
DefaultAssay(l1_teaseq_obj) <- 'ADT'
l1_wnn_clust_adt_markers <- FindAllMarkers(l1_teaseq_obj, assay = 'ADT') # using default parameters
l1_wnn_clust_adt_markers <- l1_wnn_clust_adt_markers %>%
  rename(
    l1_cluster = cluster,
    adt_epitope = gene,
    avg_log2fc_clust_vs_rest = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_clust = pct.1,
    pct_pos_in_others = pct.2
  )
l1_wnn_clust_adt_markers <- l1_wnn_clust_adt_markers %>% arrange(l1_cluster, desc(avg_log2fc_clust_vs_rest), p_val_adj)
l1_wnn_clust_adt_markers$adt_epitope <- gsub("^adt-", "", l1_wnn_clust_adt_markers$adt_epitope)
l1_wnn_clust_adt_markers <- l1_wnn_clust_adt_markers %>% relocate(c('l1_cluster','adt_epitope','avg_log2fc_clust_vs_rest','p_val_raw','p_val_adj'))
saveRDS(l1_wnn_clust_adt_markers, '/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_adt_markers.rds')

# Batch read in .rds files and write one multi-sheet XLSX spreadsheet
clust_marker_df_dir   <- '/filepath/step13_l1_teaseq_analysis/'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) |> str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
clust_marker_xlsx_sheets
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/step13_l1_teaseq_analysis/l1_wnn_clust_markers_all_assays.xlsx")  

# 19) Save L1 object and proceed to L2 subclustering of T cells ----
saveRDS(l1_teaseq_obj, '/filepath/step13_l1_teaseq_analysis/l1_teaseq_obj.rds')