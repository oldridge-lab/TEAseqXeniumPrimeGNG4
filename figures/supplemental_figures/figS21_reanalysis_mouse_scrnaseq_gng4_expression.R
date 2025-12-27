# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code for Fig. S21, panels F & G
# Fig. S21 - Sparse Gng4 RNA expression in T cells from C57BL/6J laboratory mice across diverse health and disease states.

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

# Fig S21F - LCMV-specific CD4 T cells in spleen from both acute and chronic LCMV models in C57BL/6J mice ----

# Reanalysis of Andreatta et al. A CD4+ T cell reference map delineates subtype-specific adaptation during acute and chronic viral infections. eLife. 2022. PMID 35829695. https://pubmed.ncbi.nlm.nih.gov/35829695/
# Seurat object ('ref_LCMV_CD4_mouse_release_v1.rds') obtained from Swiss Portal for Immune Cell Analysis via figshare - https://doi.org/10.6084/m9.figshare.16592693.v1

# Import Seurat object
cd4_lcmv_obj <- readRDS('/filepath/fig5/fig5_supp/ref_LCMV_CD4_mouse_release_v1.rds')

# Inspect original author annotations
table(cd4_lcmv_obj$functional.cluster) # includes 2000 Tfh_Effector and 1900 Tfh_Memory cells, as well as 6 other CD4 nonTfh clusters
Idents(cd4_lcmv_obj) <- 'functional.cluster'
DimPlot(cd4_lcmv_obj)

# Determine number of Gng4 RNA+ cells per cluster in atlas
DefaultAssay(cd4_lcmv_obj) <- 'RNA'
cd4_lcmv_obj$gng4_rna_class <- ifelse(FetchData(cd4_lcmv_obj, vars = "Gng4") > 0, "Gng4_pos", "Gng4_neg")
Idents(cd4_lcmv_obj) <- 'gng4_rna_class'
table(cd4_lcmv_obj$gng4_rna_class, cd4_lcmv_obj$functional.cluster)
# Only 20/2000 TfhEff are Gng4 RNA+ = 1%
# Even more sparse expression of Gng4 RNA in other clusters

# Find DEG between TfhEff vs TfhMem clusters
tfh_eff_vs_mem_deg <- FindMarkers(cd4_lcmv_obj, ident.1 = 'Tfh_Effector', ident.2 = 'Tfh_Memory')

# Save TfhEff vs TfhMem DEG list for Data File S8
tfh_eff_vs_mem_deg$gene <- rownames(tfh_eff_vs_mem_deg)
tfh_eff_vs_mem_deg_save <- tfh_eff_vs_mem_deg
colnames(tfh_eff_vs_mem_deg_save) <- c(
  'p_val_raw',
  'avg_log2FC_TfhEff_vs_TfhMem',
  'pct_pos_TfhEff',
  'pct_pos_TfhMem',
  'p_val_adj',
  'gene_symbol'
)
tfh_eff_vs_mem_deg_save <- tfh_eff_vs_mem_deg_save %>% arrange(p_val_adj)
write_xlsx(tfh_eff_vs_mem_deg_save, path = '/filepath/fig5/fig5_supp/mouse_tfheff_vs_tfhmem_deg.xlsx')

# Find DEG between TfhEff vs nonTfh clusters
tfh_eff_vs_nontfh <- FindMarkers(cd4_lcmv_obj, ident.1 = 'Tfh_Effector', ident.2 = c('Th1_Effector',
                                                                                     'Tcmp',
                                                                                     'Th1_Memory',
                                                                                     'Tcm',
                                                                                     'INFI_stimulated',
                                                                                     'Eomes_HI',
                                                                                     'Treg'))

# Save TfhEff vs nonTfh DEG list for Data File S8
tfh_eff_vs_nontfh$gene <- rownames(tfh_eff_vs_nontfh)
View(tfh_eff_vs_nontfh)
tfh_eff_vs_nontfh_deg_save <- tfh_eff_vs_nontfh
colnames(tfh_eff_vs_nontfh_deg_save)
colnames(tfh_eff_vs_nontfh_deg_save) <- c(
  'p_val_raw',
  'avg_log2FC_TfhEff_vs_NonTfh',
  'pct_pos_TfhEff',
  'pct_pos_NonTfh',
  'p_val_adj',
  'gene_symbol'
)
tfh_eff_vs_nontfh_deg_save <- tfh_eff_vs_nontfh_deg_save %>% arrange(p_val_adj)
write_xlsx(tfh_eff_vs_nontfh_deg_save, path = '/filepath/fig5/fig5_supp/mouse_tfheff_vs_nontfh_deg.xlsx')

# UMAP shown in Fig S21F is from SPICA visualization portal

# Fig S21G - Tumor-infiltrating lymphocyte (TIL) subsets from both B16 and MC38 tumor models in C57BL/6J mice ----

# Reanalysis of Andreatta et al. Interpretation of T cell states from single-cell transcriptomics data using reference atlases. 2021. Nat Commun. PMID 34017005. https://pubmed.ncbi.nlm.nih.gov/34017005/
# Seurat object ('ref_TILAtlas_mouse_v1.rds') obtained from Swiss Portal for Immune Cell Analysis via figshare - https://doi.org/10.6084/m9.figshare.12478571.v3

# Import Seurat object
mouse_til_obj <- readRDS('/filepath/fig5/fig5_supp/ref_TILAtlas_mouse_v1.rds')

# Inspect original author annotations
unique(mouse_til_obj$functional.cluster)
table(mouse_til_obj$functional.cluster) # 1043 Tfh
Idents(mouse_til_obj) <- 'functional.cluster' # includes Tfh cluster and 8 nonTfh CD4/CD8 T cell subsets
DimPlot(mouse_til_obj)

# Determine number of Gng4 RNA+ cells per cluster in atlas
DefaultAssay(mouse_til_obj) <- 'RNA'
mouse_til_obj$gng4_rna_class <- ifelse(FetchData(mouse_til_obj, vars = "Gng4") > 0, "Gng4_pos", "Gng4_neg")
Idents(mouse_til_obj) <- 'gng4_rna_class'
table(mouse_til_obj$gng4_rna_class, mouse_til_obj$functional.cluster)
# Only 1/1043 Tfh are Gng4 RNA+ - very sparse expression, not enriched - other clusters also very sparse

# UMAP shown in Fig S21G is from SPICA visualization portal

# No differential expression testing performed for panels A-E of Fig. S21, where Gng4 RNA expression was essentially absent in all T cell clusters