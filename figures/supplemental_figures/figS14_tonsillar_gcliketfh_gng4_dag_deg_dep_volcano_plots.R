# Sam Barnett Dubensky et al.
# Derek A. Oldridge & Laura A. Vella Labs at the Children's Hospital of Philadelphia
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned Tfh in human lymphoid tissue
# Code and data visualization for Fig. S14 (related to Fig. 5)
# Fig. S14 - GNG4+ states among the GC-like Tfh pool are distinguished by increased expression of Tfh-associated activation and gene regulatory features, related to Figure 5

# Set up working environment ----

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

# Import Seurat object from Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors)
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Define GNG4 RNA expression threshold consistent with prior figures
l3_teaseq_tfh_obj$gng4_rna_class <- ifelse(FetchData(l3_teaseq_tfh_obj, vars = "rna_GNG4") >= 1, "GNG4_RNA_high", "GNG4_RNA_low")

# Create metadata annotation column for GNG4+/GNG4- cells among tonsillar GC-like Tfh, nonGC tonsillar Tfh, and PBMC
l3_teaseq_tfh_obj$gc_gng4_class <- dplyr::case_when(
  l3_teaseq_tfh_obj$hto.tissue == "PBMC" ~ "PBMC",
  l3_teaseq_tfh_obj$hto.tissue == "Tonsil" &
    l3_teaseq_tfh_obj$gc_vs_nongc_like == "nonGC" ~ "nonGC Tonsil",
  TRUE ~ paste(
    l3_teaseq_tfh_obj$gc_vs_nongc_like,
    l3_teaseq_tfh_obj$gng4_rna_class,
    sep = "_"
  )
)
table(l3_teaseq_tfh_obj$gc_gng4_class, useNA = "ifany")

# Specify colors for GC/GNG4 annotation
gc_gng4_cols <- c(
  "GC_GNG4_RNA_high" = 'coral3',
  "GC_GNG4_RNA_low" = "#657ab0",
  "nonGC Tonsil" = "#c9a3c1",
  "PBMC" = "#e1d4c7"
)

# Order annotations for UMAP visualization
l3_teaseq_tfh_obj$gc_gng4_class <- factor(
  l3_teaseq_tfh_obj$gc_gng4_class,
  levels = c(
    "PBMC",
    "nonGC Tonsil",
    "GC_GNG4_RNA_low",
    "GC_GNG4_RNA_high"
  )
)

# Fig S14A - Schematic of GNG4 RNA+ vs RNA- tonsillar GC-like Tfh ---- 

# 3WNN UMAP
tfh_gc_gng4_class_dimplot <- DimPlot(l3_teaseq_tfh_obj, label = FALSE, reduction = "umap.wnn.harmony", order = TRUE, pt.size = 2, 
                                     group.by = 'gc_gng4_class', cols = gc_gng4_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig5/fig5_supp/tfh_gc_gng4_class_dimplot.pdf', width = 5.5, height = 4)
tfh_gc_gng4_class_dimplot
dev.off()

# Fig. S14B - DAG in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh ----

# Find DAG in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh from Signac GeneActivity analysis of ATAC modality
DefaultAssay(l3_teaseq_tfh_obj) <- 'ACT'
Idents(l3_teaseq_tfh_obj) <- 'gc_gng4_class'
gc_gng4_dag <- FindMarkers(l3_teaseq_tfh_obj, ident.1 = 'GC_GNG4_RNA_high', ident.2 = 'GC_GNG4_RNA_low', 
                           assay = 'ACT', # From GeneActivity run including large and noncoding genes
                           min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) # permissive filters to return all features for visualization
gc_gng4_dag$gene <- rownames(gc_gng4_dag)

# Set cutoff values and colors for volcano plot (P < 1e-3 & |log2FC| > 0.5)
gc_gng4_dag_volc_cols <- ifelse(
  gc_gng4_dag$avg_log2FC < -0.5 & gc_gng4_dag$p_val < 1e-3,  "#657ab0",  
  ifelse(
    gc_gng4_dag$avg_log2FC >  0.5 & gc_gng4_dag$p_val < 1e-3,  "coral3",
    "lightgrey"                            
  )
)
names(gc_gng4_dag_volc_cols)[gc_gng4_dag_volc_cols == "lightgrey"] <- "Nonsignificant"
names(gc_gng4_dag_volc_cols)[gc_gng4_dag_volc_cols == "#657ab0"] <- "Higher in GNG4-"
names(gc_gng4_dag_volc_cols)[gc_gng4_dag_volc_cols == "coral3"] <- "Higher in GNG4+"

# Find top DAG meeting thresholds to display on volcano plot
gc_gng4_dag_vol_labs_pos <- gc_gng4_dag %>% filter(avg_log2FC > 0.5) %>% filter(p_val < 1e-3) %>% arrange(p_val_adj) %>% rownames()
gc_gng4_dag_vol_labs_neg <- gc_gng4_dag %>% filter(avg_log2FC < -0.5) %>% filter(p_val < 1e-3) %>% arrange(p_val_adj) %>% rownames()

# DAG of interest meeting thresholds to label - not all will be visualized on volcano plot due to overlap
gc_gng4_dag_vol_labs_pos_list <-c(
  "PLCL2", "CXXC5", "PTPN14", "MSI2", "B3GAT1", "ASCL2", "KSR2", "SERTAD4", "TOX", "RIN3", "SMCO4", "POU2AF1", "MAFB",
  "ITGB8", "LYN", "C1QL3", "FZD3", "PDCD1", "FHL3", "IGFL2", "TNFRSF18", "PVALB","ABCA10","NCALD","NFIA-AS1","NKX2-8","CPA2",
  "CHGB", "RIMS4", "GRK3", "PDE7B", "FZD5", "KCNQ5", "GPR161", "GPR135", "TNFRSF9", "TNFRSF4", "RIN1"
)

# DAG of interest meeting thresholds to label - not all will be visualized on volcano plot due to overlap
gc_gng4_dag_vol_labs_neg_list <- c(
  "SH2D2A","MKI67", "PREX1", "NTRK1", "NUDCD1", "ENY2", "LINC00861", "PIEZO1", "SLC2A3", "HIST3H2A", "DDX21", "HIST3H2BB", "MELK",
  "HACD4", "SMTN", "MRI1", "TNNC1", "NASP", "ZNF778","ECSCR",
  "C15orf39", "KCNJ1", "HELLS", "TMEM109", "CENPF", "MAU2", "EIF4A3", "NGRN", "NUDT21", "OGFOD1"
)

# Verify DAG of interest labels are within feature list passing thresholds
setdiff(gc_gng4_dag_vol_labs_pos_list, gc_gng4_dag_vol_labs_pos)
setdiff(gc_gng4_dag_vol_labs_neg_list, gc_gng4_dag_vol_labs_neg)

# Concatenate DAG name vectors for plot
gc_gng4_dag_vol_labs <- c(gc_gng4_dag_vol_labs_pos_list, gc_gng4_dag_vol_labs_neg_list)

# Assemble volcano plot
pdf('/filepath/fig5/fig5_supp/gc_gng4_dag_volcano.pdf', height = 6, width = 9)
EnhancedVolcano(gc_gng4_dag,
                lab = rownames(gc_gng4_dag),
                selectLab = gc_gng4_dag_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'DAG IN GNG4 RNA+ VS RNA- tonsillar GC-like Tfh',
                drawConnectors = FALSE,
                pCutoff = 1e-3,
                FCcutoff = 0.5,
                pointSize = 3,
                labSize = 5,
                colCustom= gc_gng4_dag_volc_cols,
                legendPosition = 'none',
                labFace = 'italic',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  ylim(-0.1,20) +
  xlim(-1.5,2) # For visualization, restricting x-axis to omit GNG4 itself and nonsignificant features that do not meet P < 1e-3 threshold
dev.off()

# Save DAG results for GNG4+ vs GNG4- tonsillar L4 GC-like Tfh (Data S9)
gc_gng4_dag_save <- gc_gng4_dag
gc_gng4_dag_save <- gc_gng4_dag_save %>%
  rename(
    avg_log2fc_gng4_pos_vs_neg = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_gng4_pos = pct.1,
    pct_pos_in_gng4_neg = pct.2
  )
gc_gng4_dag_save <- gc_gng4_dag_save %>% arrange(p_val_adj)
gc_gng4_dag_save <- gc_gng4_dag_save %>% relocate(c('gene_symbol','avg_log2fc_gng4_pos_vs_neg','p_val_raw','p_val_adj'))
write_xlsx(gc_gng4_dag_save, '/filepath/fig5/fig5_supp/gc_gng4_dag.xlsx')

# Fig. S14C - DEG in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh ----

# Find DEG in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
Idents(l3_teaseq_tfh_obj) <- 'gc_gng4_class'
gc_gng4_deg <- FindMarkers(l3_teaseq_tfh_obj, ident.1 = 'GC_GNG4_RNA_high', ident.2 = 'GC_GNG4_RNA_low', 
                           assay = 'RNA',
                           min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) # permissive filters to return all features for visualization
gc_gng4_deg$gene <- rownames(gc_gng4_deg)

# Set cutoff values and colors for volcano plot (P < 1e-5 & |log2FC| > 0.5)
gc_gng4_deg_volc_cols <- ifelse(
  gc_gng4_deg$avg_log2FC < -0.5 & gc_gng4_deg$p_val < 1e-5,  "#657ab0",  
  ifelse(
    gc_gng4_deg$avg_log2FC >  0.5 & gc_gng4_deg$p_val < 1e-5,  "coral3",
    "lightgrey"                            
  )
)
names(gc_gng4_deg_volc_cols)[gc_gng4_deg_volc_cols == "lightgrey"] <- "Nonsignificant"
names(gc_gng4_deg_volc_cols)[gc_gng4_deg_volc_cols == "#657ab0"] <- "Higher in GNG4-"
names(gc_gng4_deg_volc_cols)[gc_gng4_deg_volc_cols == "coral3"] <- "Higher in GNG4+"

# Find top DEG meeting thresholds to display on volcano plot
gc_gng4_deg_vol_labs_pos <- gc_gng4_deg %>% filter(avg_log2FC > 0.5) %>% filter(p_val < 1e-5) %>% arrange(p_val_adj) %>% rownames()
gc_gng4_deg_vol_labs_neg <- gc_gng4_deg %>% filter(avg_log2FC < -0.5) %>% filter(p_val < 1e-5) %>% arrange(p_val_adj) %>% rownames()

# DEG of interest meeting thresholds to label - not all will be visualized on volcano plot due to overlap
gc_gng4_deg_vol_labs_pos_list <- c(
  "GNG4", "PTPN14", "KSR2", "MSI2", "DRAIC", "XXYLT1", "DAB1", "FZD3", "IKZF2", "POU2AF1", "BCL2", "NRP1", "SMCO4","P2RY8",
  "MYO6", "SGPP2", "DLEU2", "CNIH3", "RIN3", "BACH1", "B3GAT1", "WNK2", "TIGIT","SNTB1",
  "ICA1", "HIVEP1", "MYB", "PLCL2", "CXXC5", "LYN","GRIK4","LINC02099", "AC025434.1", "CAV1", "LONRF2", "DAB1-AS1", "AC021818.1", "CLDN16", "CNKSR3"
)

# DEG of interest meeting thresholds to label - not all will be visualized on volcano plot due to overlap
gc_gng4_deg_vol_labs_neg_list <- c(
  "AAK1", "RNF213", "PREX1", "KLF6", "CCND2", "CDHR3", "SMCHD1", "SLC2A3", "MAN1C1", "TSC22D3", "ATP2B4", "LINC00861", "FYB1",
  "CYTIP", "ANK3", "PABPC1", "GIMAP5", "THEMIS", "MGAT5", "ST6GALNAC3", "CD2", "ZFP36", "LINC02320", "TXK",
  "ZC3HAV1", "SAMSN1", "ACAP1", "TNIP3", "SLFN12L", "GBP5", "SARAF", "SLFN5"
)

# Verify DEG of interest labels are within feature list passing thresholds
setdiff(gc_gng4_deg_vol_labs_pos_list, gc_gng4_deg_vol_labs_pos)
setdiff(gc_gng4_deg_vol_labs_neg_list, gc_gng4_deg_vol_labs_neg)

# Concatenate DEG name label vectors for plot
gc_gng4_deg_vol_labs <- c(gc_gng4_deg_vol_labs_pos_list, gc_gng4_deg_vol_labs_neg_list)

# Assemble volcano plot
pdf('/filepath/fig5/fig5_supp/gc_gng4_deg_volcano.pdf', height = 6, width = 9)
EnhancedVolcano(gc_gng4_deg,
                lab = rownames(gc_gng4_deg),
                selectLab = gc_gng4_deg_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'DEG in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh',
                drawConnectors = FALSE,
                pCutoff = 1e-5,
                FCcutoff = 0.5,
                pointSize = 3,
                labSize = 5,
                colCustom= gc_gng4_deg_volc_cols,
                legendPosition = 'none',
                labFace = 'italic',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  xlim(-2.5,4) # For visualization, restricting x-axis to omit GNG4 itself and nonsignificant features that do not meet P < 1e-5 threshold
dev.off()

# Save DEG results for GNG4+ vs GNG4- tonsillar L4 GC-like Tfh (Data S9)
gc_gng4_deg_save <- gc_gng4_deg
gc_gng4_deg_save <- gc_gng4_deg_save %>%
  rename(
    avg_log2fc_gng4_pos_vs_neg = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_gng4_pos = pct.1,
    pct_pos_in_gng4_neg = pct.2
  )
gc_gng4_deg_save <- gc_gng4_deg_save %>% arrange(p_val_adj)
gc_gng4_deg_save <- gc_gng4_deg_save %>% relocate(c('gene_symbol','avg_log2fc_gng4_pos_vs_neg','p_val_raw','p_val_adj'))
write_xlsx(gc_gng4_deg_save, '/filepath/fig5/fig5_supp/gc_gng4_deg.xlsx')

# Fig. S14D - DEP in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh ----

# Find differentially expressed surface epitopes (ADT) in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh
DefaultAssay(l3_teaseq_tfh_obj) <- 'ADT'
Idents(l3_teaseq_tfh_obj) <- 'gc_gng4_class'
gc_gng4_dep <- FindMarkers(l3_teaseq_tfh_obj, ident.1 = 'GC_GNG4_RNA_high', ident.2 = 'GC_GNG4_RNA_low', 
                           assay = 'ADT', 
                           min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) # permissive filter to return all genes for visualization

# Remove ADT prefix from ADT column
rownames(gc_gng4_dep) <- gsub(
  "^adt-",
  "",
  rownames(gc_gng4_dep)
)
gc_gng4_dep$adt <- rownames(gc_gng4_dep)

# Set cutoff values and colors for volcano plot (P < 1e-5 & |log2FC| > 0.1)
gc_gng4_dep_volc_cols <- ifelse(
  gc_gng4_dep$avg_log2FC < -0.1 & gc_gng4_dep$p_val < 1e-5,  "#657ab0",  
  ifelse(
    gc_gng4_dep$avg_log2FC >  0.1 & gc_gng4_dep$p_val < 1e-5,  "coral3",
    "lightgrey"                            
  )
)
names(gc_gng4_dep_volc_cols)[gc_gng4_dep_volc_cols == "lightgrey"] <- "Nonsignificant"
names(gc_gng4_dep_volc_cols)[gc_gng4_dep_volc_cols == "#657ab0"] <- "Higher in GNG4-"
names(gc_gng4_dep_volc_cols)[gc_gng4_dep_volc_cols == "coral3"] <- "Higher in GNG4+"

# Find top DEP meeting thresholds to display on volcano plot
gc_gng4_dep_vol_labs_pos <- gc_gng4_dep %>% filter(avg_log2FC > 0.1) %>% filter(p_val < 1e-5) %>% arrange(p_val_adj) %>% rownames()
gc_gng4_dep_vol_labs_neg <- gc_gng4_dep %>% filter(avg_log2FC < -0.1) %>% filter(p_val < 1e-5) %>% arrange(p_val_adj) %>% rownames()

# Select DEP to visualize
gc_gng4_dep_vol_labs_pos_list <- c('CD200','CD57','CD279','TIGIT','CD304','CD172a','CD69','CD58')
# Concatenate label vectors for volcano plot
gc_gng4_dep_vol_labs <- c(gc_gng4_dep_vol_labs_pos_list, gc_gng4_dep_vol_labs_neg)

# Assemble volcano plot
pdf('/filepath/fig5/fig5_supp/gc_gng4_dep_volcano.pdf', height = 6, width = 9)
EnhancedVolcano(gc_gng4_dep,
                lab = rownames(gc_gng4_dep),
                selectLab = gc_gng4_dep_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'DEP in GNG4 RNA+ vs RNA- tonsillar GC-like Tfh',
                drawConnectors = FALSE,
                pCutoff = 1e-5,
                FCcutoff = 0.1,
                pointSize = 3,
                labSize = 5,
                colCustom= gc_gng4_dep_volc_cols,
                legendPosition = 'none',
                #labFace = 'italic',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = Inf,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) +
  xlim(-0.3,1) # For visualization, restricting x-axis to omit nonsignificant features that do not meet P < 1e-5 threshold
dev.off()

# Save DEP results for GNG4+ vs GNG4- tonsillar L4 GC-like Tfh (Data S9)
gc_gng4_dep_save <- gc_gng4_dep
gc_gng4_dep_save <- gc_gng4_dep_save %>%
  rename(
    avg_log2fc_gng4_pos_vs_neg = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_gng4_pos = pct.1,
    pct_pos_in_gng4_neg = pct.2
  )
gc_gng4_dep_save <- gc_gng4_dep_save %>% arrange(p_val_adj)
gc_gng4_dep_save <- gc_gng4_dep_save %>% relocate(c('adt','avg_log2fc_gng4_pos_vs_neg','p_val_raw','p_val_adj'))
write_xlsx(gc_gng4_dep_save, '/filepath/fig5/fig5_supp/gc_gng4_dep.xlsx')