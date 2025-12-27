# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. 3
# Fig. 3 - Multimodal GNG4 expression distinguishes activated GC-like Tfh cell states.

# Set up R working environment ----

# Set working directory to Fig 3
setwd('/filepath/fig3/')

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
library(ggrepel) # 0.9.5

# Import Seurat objects for L1-L3 TEAseq objects ----

# Import Seurat object from TEAseq Data Preprocessing Step 13 (trimodal dimensionality reduction and L1 3WNN clustering of all tonsil and peripheral blood mononuclear cells, including Harmony integration across donors)
l1_teaseq_obj <- readRDS('/filepath/step13_bulk_harmony/l1_teaseq_obj.rds')

# Import Seurat object from TEAseq Data Preprocessing Step 14 (trimodal dimensionality reduction and L2 3WNN subclustering of T cells from L1 object, including Harmony integration across donors)
l2_teaseq_tcell_obj <- readRDS('/filepath/step14_tcell_subcluster/l2_teaseq_tcell_obj.rds')

# Import Seurat object from Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors)
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Set cluster annotations for each object ----

# L1 3WNN primary cluster names
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
l1_teaseq_obj$l1_wnn_annot <- l1_teaseq_obj$wnn.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_wnn_clust_names)
l1_teaseq_obj$l1_wnn_annot <- Idents(l1_teaseq_obj)

# L2 3WNN T cell subcluster names
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
l2_teaseq_tcell_obj$l2_tcell_wnn_annot <- l2_teaseq_tcell_obj$wnn.tcell.subcluster.harmony
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'
l2_teaseq_tcell_obj <- RenameIdents(l2_teaseq_tcell_obj, l2_tcell_clust_names)
l2_teaseq_tcell_obj$l2_tcell_wnn_annot <- Idents(l2_teaseq_tcell_obj)
Idents(l2_teaseq_tcell_obj) <- 'l2_tcell_wnn_annot'

# L3 3WNN Tfh-like cell subcluster names
tfh_subclust_names <- c(
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
l3_teaseq_tfh_obj$tfh_wnn_annot <- l3_teaseq_tfh_obj$wnn.tfh.subcluster.harmony
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
l3_teaseq_tfh_obj <- RenameIdents(l3_teaseq_tfh_obj, tfh_subclust_names)
l3_teaseq_tfh_obj$tfh_wnn_annot <- Idents(l3_teaseq_tfh_obj)

# Annotate GC vs nonGC-like groups for comparisons below
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

# Specify cluster colors for visualization  ----

# L1 TEAseq object colors (primary clusters)
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

# L2 TEAseq object colors (T cell subclusters)
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

# L3 TEAseq object colors (Tfh-like subclusters)
tfh_clust_cols <- c(
  "Tfh-Circ" = '#657ab0',
  "Tcm" = "#728782",
  "Tfh-Int" = "#c49980",
  "Tfh-NFATC1" = "#c9a3c1",
  "Tfh-IL10" = "#dfcc78",
  "Tfh-Resting" = "#b6b7ba",
  "Tfh-CXCL13" = "#bad4f4",
  "Tfh-BOB1" = "#c9744d",
  "Tfh-AP1" = "#e1d4c7"
)

# L3 TEAseq GC vs nonGC-like Tfh groups
tfh_gc_vs_nongc_cols <- c(
  "Tfh-Circ" = '#2C7C7C',
  "Tcm" = "#e1d4c7", # N/A, not Tfh, specified as beige in plot and filtered in L4 Tfh-only analyses)
  "Tfh-Int" = "#E69F00",
  "Tfh-NFATC1" = "#E69F00",
  "Tfh-IL10" = "#E69F00",
  "Tfh-Resting" = "#2C7C7C",
  "Tfh-CXCL13" = "#E69F00",
  "Tfh-BOB1" = "#E69F04",
  "Tfh-AP1" = "#2C7C7C"
)

# Tissue colors
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# Fig 3A - GC vs nonGC Tfh comparison schematic ---- 

# 3WNN UMAP of L3 Tfh-like subclusters, colored by GC vs nonGC-like Tfh annotation with nonTfh Tcm in beige
tfh_gc_vs_nongc_dimplot <- DimPlot(l3_teaseq_tfh_obj, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 5.25, reduction = "umap.wnn.harmony", pt.size = 2, 
                                   group.by = 'tfh_wnn_annot', cols = tfh_gc_vs_nongc_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
pdf('/filepath/fig3/fig3a/tfh_gc_vs_nongc_dimplot.pdf', width = 5.5, height = 4)
tfh_gc_vs_nongc_dimplot
dev.off()

# Icons exported from BioRender as PDF, schematic assembled in Illustrator

# For GC vs nonGC-like Tfh comparisons in Fig 3B-C, first filter Tcm from L3 object to create Tfh-only L4 object
l4_teaseq_tfh_no_tcm_obj <- subset(l3_teaseq_tfh_obj, tfh_wnn_annot != 'Tcm')

# Fig 3B - GC vs nonGC-like Tfh group ATAC DAG Volcano Plot ---- 

# Find GC vs nonGC Tfh group differentially accessible genes (DAG) from Signac ATAC GeneActivity analysis
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ACT'
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'gc_vs_nongc_like'
gc_vs_nongc_tfh_dag <- FindMarkers(l4_teaseq_tfh_no_tcm_obj, ident.1 = 'GC', ident.2 = 'nonGC', 
                                   assay = 'ACT', # From GeneActivity run including large and noncoding genes (see parameters in GeneActivity section of TEAseq Data Preprocessing Step 15)
                                   min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, min.cells.group = 0)
gc_vs_nongc_tfh_dag$gene <- rownames(gc_vs_nongc_tfh_dag)

# Set cutoff values and colors for volcano plot (P < 0.0000000001 (1e-10) & |log2FC| > 0.5)
gc_vs_nongc_tfh_dag_volc_cols <- ifelse(
  gc_vs_nongc_tfh_dag$avg_log2FC < -0.5 & gc_vs_nongc_tfh_dag$p_val < 1e-10,  "#2C7C7C",  
  ifelse(
    gc_vs_nongc_tfh_dag$avg_log2FC >  0.5 & gc_vs_nongc_tfh_dag$p_val < 1e-10,  "#E69F04",
    "lightgrey"                            
  )
)
names(gc_vs_nongc_tfh_dag_volc_cols)[gc_vs_nongc_tfh_dag_volc_cols == "lightgrey"] <- "Nonsignificant"
names(gc_vs_nongc_tfh_dag_volc_cols)[gc_vs_nongc_tfh_dag_volc_cols == "#2C7C7C"]    <- "Higher in nonGC-like Tfh"
names(gc_vs_nongc_tfh_dag_volc_cols)[gc_vs_nongc_tfh_dag_volc_cols == "#E69F04"]    <- "Higher in GC-like Tfh"

# Find top DAG meeting thresholds to display on volcano plot (900 GC versus 42 nonGC Tfh group-enriched features)
gc_vs_nongc_tfh_dag_vol_labs_pos <- gc_vs_nongc_tfh_dag %>% filter(avg_log2FC > 0.5) %>% filter(p_val < 1e-10) %>% arrange(p_val_adj) %>% rownames()
gc_vs_nongc_tfh_dag_vol_labs_neg <- gc_vs_nongc_tfh_dag %>% filter(avg_log2FC < -0.5) %>% filter(p_val < 1e-10) %>% arrange(p_val_adj) %>% rownames()

# GC-like Tfh DAG of interest meeting thresholds to label - not all will be visualized on volcano plot due to overlap
gc_vs_nongc_tfh_dag_vol_labs_pos_list <- c(
  'AC008703.1', 'SMCO4','ASCL2','PTPN14','GNG4',
  'NCALD','MAFB','TOX','BCAT1','BTNL8',
  'PDCD1','WNK2','ADGRG1','IL21','IGFL2',
  'CEBPA','LINC01333','TSPAN5','VOPP1',
  'WSB2','BLK','POU2AF1','CXXC5','HAL',
  'THADA','SMAD1','PLCL2','CTLA4','CEP128',
  'SH2D1A','INHBB','RIN3','TSHR',
  'TOX2','ADORA2B','ADGRD1','PTPN13','CDK6',
  'CXCL13','LYN','SGPP2','ICA1','IL6R',
  'TRPS1','PIK3R3','PTPRJ','ZEB2','B3GAT1',
  'IKZF2','CNKSR3','CD58','TNFRSF18','IKZF3',
  'KSR1','KCNK5','IL6ST','FOXB1','SMAD3',
  'INSM1','LHX5','ASB1','ACOXL',
  'BATF','NR5A2','BCL2','NFATC2','SOX4',
  'ITGB8','CNIH3','IRAK3','TNC','C1QL3',
  'CD40LG','IL10','TNFSF14','ICOS',
  'RAB27A','RAB30','RAB37'
) 

# nonGC-like Tfh DAG of interest meeting thresholds to label - not all will be visualized on volcano plot due to overlap
gc_vs_nongc_tfh_dag_vol_labs_neg_list <- c(
  'KLF2','MYADM','PLAUR','CD55','SKI','GPR132','SH2D2A','DPP4','GPR157','PRKCG','SNAPC2','ADIRF'
) 

# Verify DAG of interest labels are within feature list passing significance and fold-change thresholds
gc_vs_nongc_tfh_dag_vol_labs_pos_list %in% gc_vs_nongc_tfh_dag_vol_labs_pos # All TRUE
gc_vs_nongc_tfh_dag_vol_labs_neg_list %in% gc_vs_nongc_tfh_dag_vol_labs_neg # All TRUE

# Concatenate GC and nonGC DAG vectors for volcano plot labels (using selected positive and negative fold-change genes of interest)
gc_vs_nongc_tfh_dag_vol_labs <- c(gc_vs_nongc_tfh_dag_vol_labs_pos_list, gc_vs_nongc_tfh_dag_vol_labs_neg_list)

# Assemble GC vs nonGC-like Tfh group-enriched DAG volcano plot
pdf('/filepath/fig3/fig3b/gc_vs_nongc_tfh_dag_volcano_no_tcm.pdf', height = 6.5, width = 6.5)
EnhancedVolcano(gc_vs_nongc_tfh_dag,
                lab = rownames(gc_vs_nongc_tfh_dag),
                selectLab = gc_vs_nongc_tfh_dag_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Tfh vs nonTfh-like DAG',
                drawConnectors = FALSE,
                pCutoff = 1e-10,
                FCcutoff = 0.5,
                pointSize = 2,
                labSize = 5,
                colCustom= gc_vs_nongc_tfh_dag_volc_cols,
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
  xlim(-2.1,4.15) # for visualization purposes, omit nonsignificant features that do not meet P < 1e-10 threshold
dev.off()

# Save L4 GC vs nonGC Tfh DAG results
gc_vs_nongc_tfh_dag_save <- gc_vs_nongc_tfh_dag
gc_vs_nongc_tfh_dag_save <- gc_vs_nongc_tfh_dag_save %>%
  rename(
    avg_log2fc_l4_gctfh_vs_nongctfh = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_l4gctfh = pct.1,
    pct_pos_in_l4nongctfh = pct.2
  )
gc_vs_nongc_tfh_dag_save <- gc_vs_nongc_tfh_dag_save %>% arrange(desc(avg_log2fc_l4_gctfh_vs_nongctfh), p_val_adj)
gc_vs_nongc_tfh_dag_save <- gc_vs_nongc_tfh_dag_save %>% relocate(c('gene_symbol','avg_log2fc_l4_gctfh_vs_nongctfh','p_val_raw','p_val_adj'))
saveRDS(gc_vs_nongc_tfh_dag_save, '/filepath/fig3/fig3_gcnonctfh_dag_deg_meta_supp_tables/gc_vs_nongc_tfh_dag.rds')

# Fig 3C - GC vs nonGC-like Tfh group RNA DEG Volcano Plot ---- 

# Find GC vs nonGC-like Tfh group differentially expressed genes (DEG)
l4_teaseq_tfh_no_tcm_obj <- JoinLayers(l4_teaseq_tfh_no_tcm_obj, assay = 'RNA') # join RNA layers to compute DEG
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'RNA'
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'gc_vs_nongc_like'
gc_vs_nongc_tfh_deg <- FindMarkers(l4_teaseq_tfh_no_tcm_obj, ident.1 = 'GC', ident.2 = 'nonGC', assay = 'RNA', min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, min.cells.group = 0)
gc_vs_nongc_tfh_deg$gene <- rownames(gc_vs_nongc_tfh_deg)

# Set cutoff values and colors for volcano plot (P < 0.00000000000000000001 (1e-20) & |log2FC| > 1)
gc_vs_nongc_tfh_deg_volc_cols <- ifelse(
  gc_vs_nongc_tfh_deg$avg_log2FC < -1 & gc_vs_nongc_tfh_deg$p_val < 1e-20,  "#2C7C7C",  
  ifelse(
    gc_vs_nongc_tfh_deg$avg_log2FC >  1 & gc_vs_nongc_tfh_deg$p_val < 1e-20,  "#E69F04",
    "lightgrey"                            
  )
)
names(gc_vs_nongc_tfh_deg_volc_cols)[gc_vs_nongc_tfh_deg_volc_cols == "lightgrey"] <- "Nonsignificant"
names(gc_vs_nongc_tfh_deg_volc_cols)[gc_vs_nongc_tfh_deg_volc_cols == "#2C7C7C"]    <- "Higher in nonGC-like Tfh"
names(gc_vs_nongc_tfh_deg_volc_cols)[gc_vs_nongc_tfh_deg_volc_cols == "#E69F04"]    <- "Higher in GC-like Tfh"

# Find top DEG to label on volcano plot (162 DEG in GC-like group, 76 in nonGC-like group)
gc_vs_nongc_tfh_deg_vol_labs_pos <- gc_vs_nongc_tfh_deg %>% filter(avg_log2FC > 1) %>% filter(p_val < 1e-20) %>% arrange(p_val_adj) %>% rownames()
gc_vs_nongc_tfh_deg_vol_labs_neg <- gc_vs_nongc_tfh_deg %>% filter(avg_log2FC < -1) %>% filter(p_val < 1e-20) %>% arrange(p_val_adj) %>% rownames()

# GC-like Tfh DEG of interest to label on volcano plot - some will be omitted due to overlap
gc_vs_nongc_tfh_deg_vol_labs_pos_list <- c(
  'GNG4','DRAIC','KSR2','TOX2','XXYLT1',
  'PDCD1','PTPN14','DAB1','TOX','ICA1',
  'CDK6','POU2AF1','TMEM178B','DTHD1',
  'WNK2','CEP128','IKZF2','RAB27A','SMCO4',
  'FRMD4A','NFIA','PTPRJ','P3H2','NTRK3',
  'HAL','LINC02099','MAF','MYO6','STK39',
  'IKZF3','MSI2','PTPN11','GFOD1','NPIPB4',
  'CACNA1C','PARD3','BCL6','TIGIT','RIN3',
  'MYB','EVI5','CNIH3','B3GAT1','PIK3R3',
  'CXXC5','RGS3','MYO7A','CD38','TRIB1',
  'LYN','SMAD1','CTLA4','CD200','KCNK5',
  'TEAD1','SGPP2','CD58','TNFRSF1B',
  'PTPN7'
)

# nonGC-like Tfh DEG of interest to label on volcano plot - some will be omitted due to overlap
gc_vs_nongc_tfh_deg_vol_labs_neg_list <- c(
  'GPR183','ANK3','ANXA1','IL7R','CDC14A','RORA','ZFP36L2',
  'CD55','CAMK4','PDE3B',
  'MYADM','AAK1','SCML4','KLF2','SMCHD1',
  'FOXP1','PCAT1','CCR7','VIM'
) # Other DEG of interest 'PIK3R5','PAG1','SKI','CD96','FTH1','DOCK10','RCAN3','SARAF','MYADM','RASA3','GNAQ','PLAC8','FMN1'

# Verify DEG of interest labels are within feature list passing thresholds
gc_vs_nongc_tfh_deg_vol_labs_pos_list %in% gc_vs_nongc_tfh_deg_vol_labs_pos # All TRUE
gc_vs_nongc_tfh_deg_vol_labs_neg_list %in% gc_vs_nongc_tfh_deg_vol_labs_neg # All TRUE

# Concatenate GC and nonGC DEG of interest to annotate on volcano plot
gc_vs_nongc_tfh_deg_vol_labs <- c(gc_vs_nongc_tfh_deg_vol_labs_pos_list, gc_vs_nongc_tfh_deg_vol_labs_neg_list)

# Assemble GC vs nonGC-like Tfh DEG volcano plot
pdf('/filepath/fig3/fig3c/gc_vs_nongc_tfh_deg_volcano_no_tcm.pdf', height = 6.5, width = 6.5)
EnhancedVolcano(gc_vs_nongc_tfh_deg,
                lab = rownames(gc_vs_nongc_tfh_deg),
                selectLab = gc_vs_nongc_tfh_deg_vol_labs,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'Tfh vs nonTfh-like DEG',
                drawConnectors = FALSE,
                pCutoff = 1e-20,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 5,
                colCustom= gc_vs_nongc_tfh_deg_volc_cols,
                legendPosition = 'none',
                labFace = 'italic',
                caption = NULL,
                boxedLabels = FALSE,
                parseLabels = FALSE,
                max.overlaps = 25,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
) + 
  theme(plot.subtitle = element_blank(),
        plot.title = element_blank(),
        legend.text = element_blank(),
        text = element_text(family = "sans"),
        axis.title = element_blank()) + 
  xlim(-5.2,4.1) # for visualization purposes, omit nonsignificant features that do not meet P < 1e-20 threshold
dev.off()

# Save GC vs nonGC Tfh DEG results
gc_vs_nongc_tfh_deg_save <- gc_vs_nongc_tfh_deg
gc_vs_nongc_tfh_deg_save <- gc_vs_nongc_tfh_deg_save %>%
  rename(
    avg_log2fc_l4_gctfh_vs_nongctfh = avg_log2FC,
    gene_symbol = gene,
    p_val_raw = p_val,
    pct_pos_in_l4gctfh = pct.1,
    pct_pos_in_l4nongctfh = pct.2
  )
gc_vs_nongc_tfh_deg_save <- gc_vs_nongc_tfh_deg_save %>% arrange(desc(avg_log2fc_l4_gctfh_vs_nongctfh), p_val_adj)
gc_vs_nongc_tfh_deg_save <- gc_vs_nongc_tfh_deg_save %>% relocate(c('gene_symbol','avg_log2fc_l4_gctfh_vs_nongctfh','p_val_raw','p_val_adj'))
saveRDS(gc_vs_nongc_tfh_deg_save, '/filepath/fig3/fig3_gcnonctfh_dag_deg_meta_supp_tables/gc_vs_nongc_tfh_deg_save.rds')

# Fig 3D - Meta-analysis of GC vs nonGC-like Tfh group DAG and DEG results ----

# Meta-analysis script prepared by Derek A. Oldridge, based on published method - Willer et al. (2010) Bioinformatics. METAL: fast and efficient meta-analysis of genomewide association scans. PMID 20616382.
#' Perform inverse-variance weighted meta-analysis across multiple studies:
#' This function takes matrices of p-values and effect size estimates (betas),
#' optionally with standard errors, and computes fixed-effects meta-analysis
#' statistics using inverse-variance weighting. The code is written in a vectorized
#' manner so as to run efficiently over many independent meta-analyses -- e.g.
#' separate meta-analyses by gene, integrating several modalities for each gene.
#'
#' @param pval_mat A numeric matrix of p-values. 
#'        - Each row represents a single meta-analysis unit (e.g., a gene, genetic variant, or other test).
#'        - Each column corresponds to a different study.
#'
#' @param beta_mat A numeric matrix of effect size estimates (betas).
#'        - Must have the same dimensions as pval_mat.
#'        - Each element gives the effect estimate for a study and analysis unit.
#'
#' @param stderr_mat (Optional) A numeric matrix of standard errors.
#'        - Must match dimensions of beta_mat and pval_mat if provided.
#'        - If NULL, standard errors will be estimated as beta / z-score.
#'
#' @return A data.frame with four columns:
#'         - p.meta: The meta-analysis p-value (two-tailed).
#'         - log10.p.meta: The log base 10-transformed p-meta. May be preferred/required to avoid underflow when p.meta is small.
#'         - beta.meta: The inverse-variance weighted meta-analysis beta.
#'         - stderr.meta: The combined standard error across all studies.

# Defining meta-analysis function
meta_analysis_inv_variance <- function(pval_mat, beta_mat, stderr_mat=NULL) {
  # Convert p-values and effect directions into z-scores
  # -qnorm(p/2) transforms two-tailed p-values into absolute z-scores
  # Multiplying by sign(beta) restores direction of effect
  z_mat <- -qnorm(pval_mat / 2) * sign(beta_mat)
  
  # If standard errors are not provided, estimate them using beta / z
  if(is.null(stderr_mat)) {
    stderr_mat <- beta_mat / z_mat
  }
  
  # Compute inverse-variance weights
  w_mat <- 1 / stderr_mat^2
  
  # Calculate the meta-analyzed standard error for each row (variant)
  stderr_meta <- sqrt(1 / rowSums(w_mat))
  
  # Calculate meta-analyzed beta using weighted average
  beta_meta <- rowSums(beta_mat * w_mat) / rowSums(w_mat)
  
  # Calculate meta-analysis z-score
  z_meta <- beta_meta / stderr_meta
  
  # Convert z-scores back to two-tailed p-values
  p_meta <- 2 * pnorm(-abs(z_meta))
  log10_p_meta <- (log(2) + pnorm(-abs(z_meta), log.p = TRUE)) / log(10)
  
  # Return a data frame with meta-analysis p-values (including log base 10-transformed), effect sizes, and standard errors
  return(data.frame(p.meta = p_meta, log10.p.meta = log10_p_meta, beta.meta = beta_meta, stderr.meta = stderr_meta))
}

# Prepare shared DAG vs DEG feature list for meta-analysis
nrow(gc_vs_nongc_tfh_dag) # 20338 unique annotated gene loci - note GeneActivity analysis was restricted to features annotated in RNA assay (refer to GeneActivity analysis section from TEAseq Step 15 Data Preprocessing)
nrow(gc_vs_nongc_tfh_deg) # 36601 unique transcripts, including coding and noncoding
genes_overlapping <- rownames(gc_vs_nongc_tfh_dag)[rownames(gc_vs_nongc_tfh_dag) %in% rownames(gc_vs_nongc_tfh_deg)]
length(genes_overlapping) # 20338 total features shared

# Filter and match ordering of feature sets on shared features only
deg_sorted <- gc_vs_nongc_tfh_deg[rownames(gc_vs_nongc_tfh_deg) %in% genes_overlapping,]
deg_sorted <- deg_sorted[order(rownames(deg_sorted)),]
dag_sorted <- gc_vs_nongc_tfh_dag[rownames(gc_vs_nongc_tfh_dag) %in% genes_overlapping,]
dag_sorted <- dag_sorted[order(rownames(dag_sorted)),]

# Prepare matrix summary of effect size estimates for DAG and DEG results
beta_mat <- data.frame(
  beta_deg = deg_sorted$avg_log2FC,
  beta_dag = dag_sorted$avg_log2FC
)
rownames(beta_mat) <- rownames(deg_sorted)
beta_mat <- as.matrix(beta_mat)

# Prepare matrix summary of p-values for DAG and DEG results
pval_mat <- data.frame(
  pval_deg = deg_sorted$p_val,
  pval_dag = dag_sorted$p_val
)
rownames(pval_mat) <- rownames(deg_sorted)
pval_mat <- as.matrix(pval_mat)

# Run meta-analysis function on input DAG and DEG data
meta_results <- meta_analysis_inv_variance(pval_mat = pval_mat, beta_mat = beta_mat)

# Inspect meta-analysis results
meta_results$gene <- rownames(meta_results)
meta_results <- meta_results[order(meta_results$log10.p.meta, decreasing = FALSE),]
head(meta_results)

# Extract top 10 multimodal features enriched in GC-like group based on lowest meta P-value, sorted by highest meta fold-change
top10_gc_feats_meta <- meta_results %>% filter(beta.meta > 0) %>% arrange(p.meta) %>% head(10) %>% arrange(desc(beta.meta))
top10_gc_feats_meta

# Save results for GC vs nonGC-like Tfh group DAG/DEG meta-analysis
l4_gc_vs_nongc_tfh_dag_deg_meta_results <- meta_results
l4_gc_vs_nongc_tfh_dag_deg_meta_results <- l4_gc_vs_nongc_tfh_dag_deg_meta_results %>%
  rename(
    meta_effect_size_gc_vs_nongc = beta.meta,
    gene_symbol = gene,
    meta_p_val = p.meta,
    meta_log10_p_val = log10.p.meta,
    meta_effect_size_std_err = stderr.meta
  )
l4_gc_vs_nongc_tfh_dag_deg_meta_results <- l4_gc_vs_nongc_tfh_dag_deg_meta_results %>% arrange(desc(meta_effect_size_gc_vs_nongc))
l4_gc_vs_nongc_tfh_dag_deg_meta_results <- l4_gc_vs_nongc_tfh_dag_deg_meta_results %>% relocate(c('gene_symbol','meta_effect_size_gc_vs_nongc','meta_effect_size_std_err','meta_p_val','meta_log10_p_val'))
saveRDS(l4_gc_vs_nongc_tfh_dag_deg_meta_results, '/filepath/fig3/fig3_gcnonctfh_dag_deg_meta_supp_tables/l4_gc_vs_nongc_tfh_dag_deg_meta_results.rds')

# Saving GC vs nonGC DAG, DEG, and Meta-Analysis results for supplemental data files ----

# Export all cluster marker files as sheets within one XLSX spreadsheet
clust_marker_df_dir   <- '/filepath/fig3/fig3_gcnonctfh_dag_deg_meta_supp_tables'
clust_marker_df_files <- list.files(clust_marker_df_dir, pattern = "\\.rds$", full.names = TRUE)
clust_marker_xlsx_sheets <- setNames(vector("list", length(clust_marker_df_files)), nm = basename(clust_marker_df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(clust_marker_df_files)) {
  df <- readRDS(clust_marker_df_files[i])
  clust_marker_xlsx_sheets[[i]] <- df
}
write_xlsx(clust_marker_xlsx_sheets, path = "/filepath/fig3/fig3_gcnonctfh_dag_deg_meta_supp_tables/l4_gc_vs_nongc_tfh_dag_deg_meta.xlsx")  

# Fig 3E - GNG4 ATAC and RNA 3WNN UMAP Feature Plots ----

# GNG4 ATAC GeneActivity on L3 UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'ACT'
pdf('/filepath/fig3/fig3e/gng4_atac_geneactivity_l3_umap.pdf', width = 6, height = 4)
FeaturePlot(l3_teaseq_tfh_obj, label = FALSE, reduction = "umap.wnn.harmony", pt.size = 2, features = 'GNG4', order = TRUE, 
            cols = c('#e2e2e2','coral2')) + coord_fixed() + NoAxes() + 
  theme(
    plot.title = element_blank(),
    legend.position      = c(0.925, 0.5),
    legend.justification = c(0, 0.5),
    legend.margin        = margin(2, 2, 2, 2)
  )
dev.off()

# GNG4 RNA on L3 UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig3/fig3e/gng4_rna_l3_umap.pdf', width = 6, height = 4)
FeaturePlot(l3_teaseq_tfh_obj, label = FALSE, reduction = "umap.wnn.harmony", pt.size = 2, features = 'GNG4', order = TRUE, 
            cols = c('#e2e2e2','coral2')) + coord_fixed() + NoAxes() + 
  theme(
    plot.title = element_blank(),
    legend.position      = c(0.925, 0.5),
    legend.justification = c(0, 0.5),
    legend.margin        = margin(2, 2, 2, 2)
  )
dev.off()

# Fig 3F - L1 TEAseq 3WNN UMAP with Tfh specificity of GNG4 RNA  ----

# GNG4 RNA expression on L1 TEAseq 3WNN UMAP
pdf('/filepath/fig3/fig3f/l1_wnn_umap_gng4_rna.pdf', width = 6, height = 6)
DefaultAssay(l1_teaseq_obj) <- 'RNA'
FeaturePlot(l1_teaseq_obj, label = FALSE, reduction = "umap.wnn.harmony", pt.size = 1, features = 'GNG4', order = TRUE, 
            cols = c('#e2e2e2','coral2')) + coord_fixed() + NoAxes() + 
  theme(
    plot.title = element_blank()
  )
dev.off()

# For reference, annotated L1 3WNN TEAseq object as in Fig 1B
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
l1_teaseq_obj$l1_wnn_annot <- l1_teaseq_obj$wnn.bulk.cluster.harmony
Idents(l1_teaseq_obj) <- 'l1_wnn_annot'
l1_teaseq_obj <- RenameIdents(l1_teaseq_obj, l1_wnn_clust_names)
l1_teaseq_obj$l1_wnn_annot <- Idents(l1_teaseq_obj)
table(l1_teaseq_obj$l1_wnn_annot)
DimPlot(l1_teaseq_obj, reduction = "umap.wnn.harmony", label = TRUE, label.box = TRUE, label.size = 6, pt.size = 0.5, repel = TRUE, group.by = 'l1_wnn_annot', cols = l1_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())

# Fig 3G - L2 TEAseq 3WNN UMAP with Tfh specificity of GNG4 RNA ----

# GNG4 RNA expression on L2 TEAseq 3WNN UMAP
pdf('/filepath/fig3/fig3g/l2_tcell_wnn_umap_gng4_rna.pdf', width = 6, height = 6)
DefaultAssay(l2_teaseq_tcell_obj) <- 'RNA'
FeaturePlot(l2_teaseq_tcell_obj, label = FALSE, reduction = "umap.wnn.harmony", pt.size = 1.25, features = 'GNG4', order = TRUE, 
            cols = c('#e2e2e2','coral2')) + coord_fixed() + NoAxes() + 
  theme(
    plot.title = element_blank()
  )
dev.off()

# For reference L2 3WNN UMAP with T cell subclusters annotated as in Fig 1H
Idents(l2_teaseq_tcell_obj) <- 'tcell_wnn_annot'
tcell_wnn_umap_clust_dimplot <- DimPlot(l2_teaseq_tcell_obj, reduction = "umap.wnn.harmony", pt.size = 0.75, group.by = 'tcell_wnn_annot', cols = tcell_clust_cols) + coord_fixed() + NoLegend() + NoAxes() + theme(plot.title = element_blank())
tcell_wnn_umap_coords <- FetchData(l2_teaseq_tcell_obj, vars = c("harmonywnnUMAP_1", "harmonywnnUMAP_2")) %>%
  mutate(tcell_wnn_annot = Idents(l2_teaseq_tcell_obj),
         fill_col = tcell_clust_cols[as.character(tcell_wnn_annot)])
l1_wnn_umap_lab_df <- tcell_wnn_umap_coords %>% 
  group_by(tcell_wnn_annot, fill_col) %>% 
  summarize(UMAP_1 = median(harmonywnnUMAP_1), UMAP_2 = median(harmonywnnUMAP_2), .groups = 'drop')
tcell_wnn_umap_clust_lab <- tcell_wnn_umap_clust_dimplot +
  new_scale_color() + 
  geom_label_repel(
    data = l1_wnn_umap_lab_df,
    aes(x = UMAP_1, y = UMAP_2, label = tcell_wnn_annot, fill = fill_col, color = 'black'),
    size = 5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines"),
    force = 1,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  NoLegend()
tcell_wnn_umap_clust_lab

# Fig 3H - Connected dot plots for mean ADT expression in GNG4 RNA+ vs RNA- Tfh groups per tonsil donor ----

# Specify threshold for GNG4 RNA+ versus GNG4 RNA- cells in L4 Tfh-only object
l4_teaseq_tfh_no_tcm_obj$gng4_rna_class <- ifelse(FetchData(l4_teaseq_tfh_no_tcm_obj, vars = "rna_GNG4") >= 1, "GNG4_RNA_high", "GNG4_RNA_low")
table(l4_teaseq_tfh_no_tcm_obj$gng4_rna_class) # 2996 total Tfh - 597 RNA+ (19.9%) vs 2399 (80.1%)
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'gng4_rna_class'

# Get normalized ADT expression values between GNG4 RNA high versus low groups per tonsil donor
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ADT'
adt_feature <- rownames(l4_teaseq_tfh_no_tcm_obj) # consider all features in ADT assay
tfh_notcm_gng4class_adt_md <- l4_teaseq_tfh_no_tcm_obj %>% # using L4 Tfh-only object (i.e. Tcm cluster excluded)
  subset(hto.tissue == "Tonsil") %>% # only tonsil samples
  FetchData(vars = c("gng4_rna_class", "hto.donor", adt_feature)) %>% # Get GNG4 RNA class, donor, and ADT feature values
  rename(
    Class = gng4_rna_class,                                      
    Donor = hto.donor
  )

# Get mean ADT signal per donor within each GNG4 RNA class
donor_averages <- tfh_notcm_gng4class_adt_md %>%
  group_by(Class, Donor) %>%
  summarise(
    across(starts_with("adt-"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    .groups = "drop"
  )

# Transpose output ADT data, import into Prism for statistical analysis via paired t-test with two-stage step-up procedure of Benjamini, Krieger, and Yekutieli
donor_averages_t <- as.data.frame(t(donor_averages))
donor_averages_t$adt <- rownames(donor_averages_t)
write_xlsx(donor_averages_t, '/filepath/fig3/fig3h/adt_avg_ton_gng4_subsets_df_transpose.xlsx')

# Fig 3I - Flow contour plots showing gating strategy of Tfh and Gy4 histogram ----

# PDF exported from FlowJo, further annotated in Illustrator

# Fig 3J -  Gy4+ percentage in Tfh subsets vs nonTfh and Naive CD4 T cells ----

# Tfh Susbets % Gy4+ Barplot
tfh_gng4_flow_df <- read_excel("/filepath/fig3/fig3_flow_data.xlsx", sheet = "Gy4_Freq")

tonsil_barplot_cols <- c(
  "Naïve" = '#d6d6d6',
  "nonTfh" = "#824100",
  "PD1-" = "#ac5600",
  "PD1+" = "#d66b00",
  "PD1br" = "#ff7f00"
)

pbmc_barplot_cols <- c(
  "Naïve" = '#d6d6d6',
  "nonTfh" = "#002c58",
  "PD1-" = "#004182",
  "PD1+" = "#0056ac",
  "PD1br" = "#007fff"
)

tfh_gng4_flow_df_long <- tfh_gng4_flow_df %>%
  pivot_longer(
    cols      = Naïve:PD1br,               
    names_to  = "Subset",                  
    values_to = "Freq"                     
  )

tfh_gng4_flow_df_long <- tfh_gng4_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c("Naïve", "nonTfh", "PD1-", "PD1+", "PD1br")
    )
  )

# Tonsil sample percentage Gy4+ barplot with 95% confidence intervals
tfh_gng4_flow_df_ton <- tfh_gng4_flow_df_long %>% subset(Tissue == 'Tonsil')
ton_tfh_gy4_freq_barplot <- ggplot(tfh_gng4_flow_df_ton, aes(x = Subset, y = Freq, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) +
  geom_jitter(
    data      = subset(tfh_gng4_flow_df_ton, Tissue == "Tonsil"),
    shape     = 21,
    color   = "black",
    stroke  = 1.5,
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.6,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_blank(),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = tonsil_barplot_cols) +
  theme(axis.title = element_blank()) +
  scale_y_continuous(breaks = seq(0, 60, 5)) + 
  coord_cartesian(ylim = c(0, 60))
ton_tfh_gy4_freq_barplot
pdf('/filepath/fig3/fig3j/ton_tfh_gy4_freq_barplot_2.pdf', width = 5, height = 5)
ton_tfh_gy4_freq_barplot
dev.off()

# PBMC sample percentage Gy4+ barplot with 95% confidence intervals
tfh_gng4_flow_df_pbmc <- tfh_gng4_flow_df_long %>% subset(Tissue == 'PBMC')
pbmc_tfh_gy4_freq_barplot <- ggplot(tfh_gng4_flow_df_pbmc, aes(x = Subset, y = Freq, fill = Subset)) +
  stat_summary(
    fun    = mean, 
    geom   = "bar",
    color  = "black",
    width  = 0.7,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "errorbar",
    width     = 0.2,
    size      = 0.5
  ) +
  geom_jitter(
    data      = subset(tfh_gng4_flow_df_pbmc, Tissue == "PBMC"),
    shape     = 21,
    color   = "black",
    stroke  = 1.5,
    width   = 0.15,
    size    = 3.5,
    alpha   = 0.6,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y      = element_blank(),
    axis.title.x = element_blank(),
    axis.title       = element_text(size = 20)) +
  scale_fill_manual(values = pbmc_barplot_cols) +
  theme(axis.title = element_blank()) +
  scale_y_continuous(breaks = seq(0, 60, 5)) +
  coord_cartesian(ylim = c(0, 60))
pbmc_tfh_gy4_freq_barplot
pdf('/filepath/fig3/fig3j/pbmc_tfh_gy4_freq_barplot.pdf', width = 5, height = 5)
pbmc_tfh_gy4_freq_barplot
dev.off()

# Data file S10 - 95% confidence interval calculations for Gy4 percentage and GMFI values in Tfh versus nonTfh subsets of interest ----

# Select sheets from raw flow data spreadsheet
xlsx_path <- "/filepath/fig3/fig3_flow_data.xlsx"
sheets <- c("Gy4_GMFI", "Gy4_Freq")

# Function to retrieve tissue and donor groups from raw data, then compute relevant 95CI metrics
process_sheet <- function(sheet_name) {
  df <- read_excel(xlsx_path, sheet = sheet_name)
  subset_cols <- df %>%
    select(-any_of(c("Tissue","Donor"))) %>%
    select(where(is.numeric)) %>%
    names()
  
  df %>%
    pivot_longer(all_of(subset_cols), names_to = "Subset", values_to = "Value") %>%
    group_by(Tissue, Subset) %>%
    summarise(
      n = n(),
      mean = mean(Value),
      sd = sd(Value),
      se = sd / sqrt(n),
      ci = qt(0.975, df = n - 1) * se,
      ymin = mean - ci,
      ymax = mean + ci,
      .groups = "drop"
    ) %>%
    mutate(sheet = sheet_name, .before = 1)
}

# Execute function
all_stats <- map_dfr(sheets, process_sheet)

# Save 95CI metrics
writexl::write_xlsx(
  list(Stats_by_Group = all_stats %>%
         select(sheet, Tissue, Subset, n, mean, sd, se, ci, ymin, ymax)),
  path = "/filepath/fig3/fig3j/fig3j_tfh_gy4_freq_gmfi_95ci_stats.xlsx"
)

all_stats_t <- as.data.frame(t(all_stats))
all_stats_t$variables <- rownames(all_stats_t)
all_stats_t <- all_stats_t %>% select(variables, everything())
writexl::write_xlsx(
  all_stats_t,
  path = "/filepath/fig3/fig3j/fig3j_tfh_gy4_freq_gmfi_95ci_stats_transpose.xlsx"
)

# Fig 3K - Flow contour plots showing features correlated versus anticorrelated with Gy4 protein expression in tonsil Tfh ----

# PDF exported from FlowJo, further annotated in Illustrator 

# Fig 3L-3N - Connected dot plots of feature expression percentage in Gy4+ vs Gy4- Tfh groups per tonsil donor ----

# Percentage values exported from FlowJo, analyzed and visualized in Prism