# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for Fig. 5
# Fig. 5. GNG4 expression in Tfh is associated with activation in vitro and in vivo.

# Set up R working environment ----

# Set working directory to Fig 5
setwd('/filepath/fig5/')

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

# Prepare L3 Tfh-like subclustering and L4 Tfh-only filtering Seurat objects with cluster annotations and colors for visualization ----

# Import Seurat object from Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors)
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Annotate L3 Tfh-like cell clusters
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

# Annotate L3 GC vs nonGC-like Tfh groups
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

# Filter Tcm from L3 object to create Tfh-only L4 object
ncol(l3_teaseq_tfh_obj) # 3657 cells total, including Tcm
l3_teaseq_tfh_obj@meta.data %>% summarise(x = sum(tfh_wnn_annot == 'Tcm')) %>% print() # 661 Tcm
l4_teaseq_tfh_no_tcm_obj <- subset(l3_teaseq_tfh_obj, tfh_wnn_annot != 'Tcm')
ncol(l4_teaseq_tfh_no_tcm_obj) # 2996 cells remaining in L4 Tfh-only object

# Assign colors to L3 clusters for visualization
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

# Assign colors to tissues
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# Fig 5A - ATAC coverage plot of GNG4 locus with DAPs across L3 clusters ---- 

# Find GNG4 DAP in GC vs nonGC-like Tfh groups (L4 Tfh subset analysis, Tcm excluded)
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ATAC'
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'gc_vs_nongc_like'
gc_vs_nongc_tfh_dap <- FindMarkers(l4_teaseq_tfh_no_tcm_obj, ident.1 = 'GC', ident.2 = 'nonGC', assay = 'ATAC', 
                                   min.pct = 0, logfc.threshold = 0) # default filters relaxed initially to explore all peaks
gc_vs_nongc_tfh_dap_gr <- StringToGRanges(rownames(gc_vs_nongc_tfh_dap), sep = c("-", "-"))
gc_vs_nongc_tfh_dap_closest_feat <- ClosestFeature( # annotating closest gene to each peak
  object = l4_teaseq_tfh_no_tcm_obj,
  regions    = gc_vs_nongc_tfh_dap_gr,
  annotation = Annotation(l4_teaseq_tfh_no_tcm_obj)
)
gc_vs_nongc_tfh_dap$closest_gene <- gc_vs_nongc_tfh_dap_closest_feat$gene_name
gc_vs_nongc_tfh_dap$transcript_id <- gc_vs_nongc_tfh_dap_closest_feat$tx_id
gc_vs_nongc_tfh_dap$atac_peak <- rownames(gc_vs_nongc_tfh_dap)

# Filter for statistically significant GNG4 DAPs in GC vs nonGC-like Tfh groups and save for data file S7
gc_vs_nongc_tfh_dap_save <- gc_vs_nongc_tfh_dap %>%
  rename(
    avg_log2FC_l4_tfh_gc_vs_nongc = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_gctfh = pct.1,
    pct_pos_in_nongctfh = pct.2
  ) %>%
  filter(p_val_adj < 0.05)
gc_vs_nongc_tfh_dap_save <- gc_vs_nongc_tfh_dap_save %>% arrange(desc(avg_log2FC_l4_tfh_gc_vs_nongc), p_val_adj)
gc_vs_nongc_tfh_dap_save <- gc_vs_nongc_tfh_dap_save %>% relocate(c('atac_peak','closest_gene','avg_log2FC_l4_tfh_gc_vs_nongc','p_val_raw','p_val_adj'))
saveRDS(gc_vs_nongc_tfh_dap_save, "/filepath/fig5/fig5a_l4_gc_vs_nongc_tfh_atac_dap_padj05filt.rds")

# Annotate GNG4 DAP of interest for coverage plot visualization in Fig 5A
# Filtered dataframe for pvaladj < 0.05 and open in > 10% of cells
# All GNG4 DAP map to transcript variant ENST00000450593, GNG4-204, 4978bp, 75aa
# Refer to Supplementary Methods for DAP region annotation - visualized DAP regions in hg38 using UCSC Genome Browser with track for ENCODE Candidate Cis-Regulatory Elements (cCREs across human cell types) 
# chr1-235649534-235650662 # Centers around two promoter-like sequence elements (one CTCF-bound), contains partial ends of pELS elements (one CTCF-bound) 
# chr1-235648709-235649522 # Centers around proximal enhancer-like sequence element (CTCF-bound), contains partial ends of two other pELS elements (one CTCF-bound)
# chr1-235641919-235642869 # Centers around distal enhancer-like sequence element (CTCF-bound), contains another complete dELS, partial end of a third dELS

# Create GRanges for GNG4 DAP of interest
gng4_dap_df <- data.frame(
  seqnames = "chr1",
  start = c(235641919, 235649534, 235648709), # start coordinates for dELS, PLS, and pELS-containing GNG4 DAP regions of interest
  end   = c(235642869, 235650662, 235649522) # end coordinates for dELS, PLS, and pELS-containing GNG4 DAP regions of interest
)
gng4_dap_gr <- makeGRangesFromDataFrame(gng4_dap_df)
rdbu_colors <- brewer.pal(11, "RdBu")
gng4_dap_gr$color <- c('darkgrey',rdbu_colors[3],'darkgrey')

# Infer links in GNG4 peak accessibility to RNA expression
DefaultAssay(l3_teaseq_tfh_obj) <- 'ATAC'
l3_teaseq_tfh_obj <- RegionStats(l3_teaseq_tfh_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
l3_teaseq_tfh_obj <- JoinLayers(l3_teaseq_tfh_obj, assay = 'RNA')
l3_teaseq_tfh_obj <- LinkPeaks(
  object = l3_teaseq_tfh_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = 'GNG4')

# Save GNG4 links for data file S7
l3_tfh_gng4_links_df <- as.data.frame(Links(l3_teaseq_tfh_obj))
l3_tfh_gng4_links_df$strand <- NULL # for ATAC +/- strand information is not defined
l3_tfh_gng4_links_df$seqnames <- NULL # chr1 for all
l3_tfh_gng4_links_df <- l3_tfh_gng4_links_df %>%
  rename(
    link_start_coord = start,
    link_end_coord = end,
    link_width = width,
    linked_atac_peak = peak,
    linked_gene = gene,
    rna_expr_pearson_corr = score,
    pearson_corr_z_score = zscore,
    z_score_p_val = pvalue
  )
l3_tfh_gng4_links_df <- l3_tfh_gng4_links_df %>% arrange(desc(rna_expr_pearson_corr))
l3_tfh_gng4_links_df <- l3_tfh_gng4_links_df %>% relocate(c('linked_gene','linked_atac_peak','link_start_coord','link_end_coord','link_width','rna_expr_pearson_corr','pearson_corr_z_score','z_score_p_val'))
View(l3_tfh_gng4_links_df)
saveRDS(l3_tfh_gng4_links_df, "/filepath/fig5/fig5a_l3_tfh_gng4_dap_rna_link_stats.rds")

# Set Tfh identity order for coverage plot
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
wnn_tfh_ident_order <- c('Tcm','Tfh-AP1','Tfh-Circ','Tfh-Resting','Tfh-IL10','Tfh-Int','Tfh-NFATC1','Tfh-CXCL13','Tfh-BOB1')
Idents(l3_teaseq_tfh_obj) <- factor(
  Idents(l3_teaseq_tfh_obj),
  levels = rev(wnn_tfh_ident_order)
)

# Setting all clusters to same color for visualization
tfh_coverage_cols <- c(
  "Tfh-Circ" = 'black',
  "Tcm" = "black",
  "Tfh-Int" = "black",
  "Tfh-NFATC1" = "black",
  "Tfh-IL10" = "black",
  "Tfh-Resting" = "black",
  "Tfh-CXCL13" = "black",
  "Tfh-BOB1" = "black",
  "Tfh-AP1" = "black"
)

# Assemble GNG4 locus coverage plot with annotation
DefaultAssay(l3_teaseq_tfh_obj) <- 'ATAC'
gng4_cov_plot <- CoveragePlot(l3_teaseq_tfh_obj, region = 'chr1-235641500-235651100', annotation = TRUE, # coordinates zoomed to three DAP of interest within GNG4 Locus
                              peaks = FALSE, 
                              links = TRUE, region.highlight = gng4_dap_gr, 
                              heights = c(7, 1)) & scale_fill_manual(values = tfh_coverage_cols)
gng4_cov_plot[[1]][[1]] <- gng4_cov_plot[[1]][[1]] + theme(axis.ticks = element_blank()) + labs(y = 'Normalized Peak Signal (0-97)') # y-axis label renamed for visualization, range first verified before renaming
gng4_cov_plot[[1]][[2]] <- gng4_cov_plot[[1]][[2]] + scale_color_manual(values = c('black','black'))
gng4_cov_plot_links <- gng4_cov_plot[[1]][[3]]
gng4_cov_plot_links <- gng4_cov_plot_links +
  scale_color_gradient(
    low  = "lightgrey",
    high = rdbu_colors[3],
    name = "RNA Linkage" 
  ) +
  theme(
    legend.key.size   = unit(0.3, "cm"),
    axis.text.x = element_text(size = 5),
    axis.title.x = element_blank(),
    legend.text       = element_text(size = 8),
    legend.title      = element_blank()
  )
gng4_cov_plot[[1]][[3]] <- gng4_cov_plot_links
pdf('/filepath/fig5/fig5a/tfh_gng4_link_plot_zoom_annotation.pdf', width = 5, height = 4)
gng4_cov_plot
dev.off()

# Create adjoining GNG4 RNA dotplot, scaled across clusters 
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
Idents(l3_teaseq_tfh_obj) <- 'tfh_wnn_annot'
wnn_tfh_ident_order <- c('Tcm','Tfh-AP1','Tfh-Circ','Tfh-Resting','Tfh-IL10','Tfh-Int','Tfh-NFATC1','Tfh-CXCL13','Tfh-BOB1')
Idents(l3_teaseq_tfh_obj) <- factor(
  Idents(l3_teaseq_tfh_obj),
  levels = wnn_tfh_ident_order
)
wnn_tfh_gng4_dotplot_feats <- c('GNG4')
wnn_tfh_gng4_dotplot <- DotPlot(l3_teaseq_tfh_obj, features = wnn_tfh_gng4_dotplot_feats,
                                   cluster.idents = FALSE, cols = 'RdBu', scale = TRUE, scale.by = 'size', dot.scale = 12) + 
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.justification.right = 'center',
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8, 'sans'),
        legend.margin = margin(t = 2, b = 2),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 0.1)) +
  guides(
    size = guide_legend(
      title = "",
      order = 1,
      reverse = TRUE,
      barwidth = 0.25
    ),
    color = guide_colorbar(
      title = "",
      order = 2,
      barwidth = 2,
      barheight = 8
    )
  )  + 
  scale_color_gradient2(
    low = rdbu_colors[11],
    mid = "white",
    high = rdbu_colors[3],
    midpoint = 0
  )
pdf('/filepath/fig5/fig5a/wnn_tfh_gng4_dotplot.pdf', height = 6, width = 2.75)
wnn_tfh_gng4_dotplot
dev.off()

# Fig 5B - Scatter plot of TF Motif Enrichment (ATAC-based) in DAPs of GNG4 RNA+ Tfh versus RNA- Tfh ----

# Define threshold for GNG4 RNA+ vs RNA- Tfh (L4 Tfh-only object, same threshold as used in Fig 3H)
l4_teaseq_tfh_no_tcm_obj$gng4_rna_class <- ifelse(FetchData(l4_teaseq_tfh_no_tcm_obj, vars = "rna_GNG4") >= 1, "GNG4_RNA_high", "GNG4_RNA_low")
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'gng4_rna_class'

# Find ATAC DAPs between GNG4 RNA+ vs RNA- Tfh
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ATAC'
gng4_rna_class_dap <- FindMarkers(l4_teaseq_tfh_no_tcm_obj, ident.1 = 'GNG4_RNA_high', ident.2 = 'GNG4_RNA_low', assay = 'ATAC') # using default logfc.threshold and min.pct values as we will filter in following step
gng4_rna_class_dap <- gng4_rna_class_dap %>% filter(p_val_adj < 0.05)
gng4_rna_class_dap_gr <- StringToGRanges(rownames(gng4_rna_class_dap), sep = c("-", "-"))
gng4_rna_class_dap_closest_feat <- ClosestFeature( # annotating DAPs by finding closest annotated gene
  object = l4_teaseq_tfh_no_tcm_obj,
  regions    = gng4_rna_class_dap_gr,
  annotation = Annotation(l4_teaseq_tfh_no_tcm_obj)
)
gng4_rna_class_dap$closest_gene <- gng4_rna_class_dap_closest_feat$gene_name
gng4_rna_class_dap$atac_peak <- rownames(gng4_rna_class_dap)
View(gng4_rna_class_dap)

# Save ATAC DAPs between GNG4 RNA+ vs RNA- L4 Tfh for Data File S7
gng4_rna_class_dap_save <- gng4_rna_class_dap %>%
  rename(
    avg_log2FC_l4tfh_gng4rna_pos_vs_neg = avg_log2FC,
    p_val_raw = p_val,
    pct_pos_in_gng4rna_pos_tfh = pct.1,
    pct_pos_in_gng4rna_neg_tfh = pct.2
  ) %>%
  filter(p_val_adj < 0.05)
gng4_rna_class_dap_save <- gng4_rna_class_dap_save %>% arrange(desc(avg_log2FC_l4tfh_gng4rna_pos_vs_neg), p_val_adj)
gng4_rna_class_dap_save <- gng4_rna_class_dap_save %>% relocate(c('atac_peak','closest_gene','avg_log2FC_l4tfh_gng4rna_pos_vs_neg','p_val_raw','p_val_adj'))
saveRDS(gng4_rna_class_dap_save, "/filepath/fig5/fig5b_l4_tfh_gng4_rna_pos_vs_neg_atac_daps_pval05_filt.rds")

# Before performing TF motif enrichment analysis, filter for DAP regions that are 1) enriched in GNG4 RNA+ Tfh rather than RNA-, and 2) statistically significant by Bonferroni-adjusted Wilcoxon P < 0.05
gng4_rna_class_dap <- gng4_rna_class_dap %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05) 
gng4_rna_class_dap_list <- rownames(gng4_rna_class_dap)
  
# Determine background peak set for motif enrichment analysis
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ATAC'
tfh_all_open_peaks <- AccessiblePeaks(l4_teaseq_tfh_no_tcm_obj) # get all open peaks in L4 Tfh-only object for peak matching to the GNG4 RNA+ Tfh DAP set
tfh_atac_meta_feat <- GetAssayData(l4_teaseq_tfh_no_tcm_obj, assay = "ATAC", layer = "meta.features")
gng4_rna_class_dap_bg_match <- MatchRegionStats( # finding background peak set with matching technical variables such as G/C nucleotide content
  meta.feature = tfh_atac_meta_feat[tfh_all_open_peaks, ],
  query.feature = tfh_atac_meta_feat[gng4_rna_class_dap_list, ]
)

# Perform TF motif enrichment on DAP regions enriched in GNG4 RNA+ Tfh group
gng4_rna_class_dap_motif_enr <- FindMotifs(
  object = l4_teaseq_tfh_no_tcm_obj,
  features = gng4_rna_class_dap_list,
  background = gng4_rna_class_dap_bg_match # background peak set determined above
)

# Determine -log10 adjusted p-values for scatter plot visualization
gng4_rna_class_dap_motif_enr$log10padj <- -log10(gng4_rna_class_dap_motif_enr$p.adjust)

# Set significance and fold change cutoffs for scatter plot
gng4_rna_class_dap_motif_enr$pass_threshold <- ifelse(
  gng4_rna_class_dap_motif_enr$p.adjust < 1e-10 & gng4_rna_class_dap_motif_enr$fold.enrichment > 1.25,
 'yes',
 'no'
)

# For Fig 5B visualization, only labeling motifs that passed thresholds, and further raising fold-enrichment threshold for visualization
gng4_rna_class_dap_motif_enr$plot_label <- ifelse(
  gng4_rna_class_dap_motif_enr$log10padj > 10 & 
    gng4_rna_class_dap_motif_enr$fold.enrichment > 1.33, # threshold is 1.25 but we label only > 1.3 here for visualization
  gng4_rna_class_dap_motif_enr$motif.name,
  NA
)

# Finding number of TF motifs passing defined enrichment thresholds (indicated in Fig 5C schematic)
gng4_rna_class_dap_motif_enr %>% filter(p.adjust < 1e-10) %>% filter(fold.enrichment > 1.25) %>% 
  rownames() %>% length() %>% print() # 63 total enriched TF passing thresholds
View(gng4_rna_class_dap_motif_enr)

# Save L4 Tfh GNG4 RNA+ DAP region TF motif enrichment results for Data File S7
gng4_rna_class_dap_motif_enr_save <- gng4_rna_class_dap_motif_enr %>%
  rename(
    tf_motif_jaspar2020 = motif,
    tf_motif_name = motif.name
  )
gng4_rna_class_dap_motif_enr_save$plot_label <- NULL
gng4_rna_class_dap_motif_enr_save <- gng4_rna_class_dap_motif_enr_save %>% arrange(desc(fold.enrichment),p.adjust)
gng4_rna_class_dap_motif_enr_save <- gng4_rna_class_dap_motif_enr_save %>% relocate(c('tf_motif_jaspar2020','tf_motif_name','fold.enrichment'))
View(gng4_rna_class_dap_motif_enr_save)
saveRDS(gng4_rna_class_dap_motif_enr_save, "/filepath/fig5/fig5b_l4_tfh_gng4_rna_pos_vs_neg_atac_daps_motif_enrichment.rds")

# Define flag to label enriched motifs in scatter plot
gng4_rna_class_dap_motif_enr$pass_threshold <- with(
  gng4_rna_class_dap_motif_enr,
  log10padj > 10 & fold.enrichment > 1.25
)

# Assemble motif enrichment scatter plot
pdf('/filepath/fig5/fig5b/tfh_dap_motif_enr_redcol_signif.pdf', width = 8.5, height = 6)
ggplot(gng4_rna_class_dap_motif_enr, aes(x = fold.enrichment, y = log10padj)) +
  geom_point( 
    data  = subset(gng4_rna_class_dap_motif_enr, !pass_threshold), # motifs below thresholds colored light grey
    shape = 21, size = 2.5, stroke = 0.5, alpha = 0.75,
    color = "black", fill = "lightgrey",
    show.legend = FALSE
  ) +
  geom_point(
    data  = subset(gng4_rna_class_dap_motif_enr, pass_threshold),
    shape = 21, size = 2.5, stroke = 0.5,
    color = "black", fill = '#D6604D', # motifs above thresholds colored coral
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(gng4_rna_class_dap_motif_enr, pass_threshold),  # only names for motifs above thresholds are labeled
    aes(label = plot_label),
    size = 4, max.overlaps = 50, show.legend = FALSE
  ) +
  labs(
    y = expression(-log[10]~adj.~italic(p)),
    x = 'Motif Fold Enrichment / Matched Background Peaks',
    title = NULL
  ) +
  theme_bw(base_size = 14, base_family = 'sans') +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14)
  ) +
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'black', linewidth = 0.5) + # -log10 adjusted P-value threshold
  geom_vline(xintercept = 1.25, linetype = 'dashed', color = 'black', linewidth = 0.5) # fold-enrichment threshold
dev.off()

# Fig 5C - Inferring regulatory TF for GNG4 RNA expression in Tfh from SCENIC and ATAC-based motif enrichment approaches ----

# Import SCENIC analysis results for L3 object after complete pipeline run
setwd('/filepath/l3_tfh_scenic_output_folder/')
scenicOptions <- readRDS("/filepath/l3_tfh_scenic_output_folder/scenicOptions4.rds") # scenicOptions object
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo") # list of TF and top-scoring motifs in v10 database for predicted target genes

# Filter regulonTargetsInfo output for TF regulons including GNG4 (columns renamed in Data File S7 for clarity)
gng4_regulons <- regulonTargetsInfo %>%
  filter(gene == 'GNG4')
saveRDS(gng4_regulons, "/filepath/fig5/fig5c_l3_tfh_scenic_regulonTargetsInfo_output_gng4.rds")

# Schematic assembled in Illustrator shows number of TFs inferred to regulate GNG4 RNA expression by SCENIC (RNA-based) versus TF motif enrichment approaches (ATAC-based)
# 13 TF regulons from SCENIC analysis - full script provided in SCENIC folder corresponding to L3 Tfh-like subclustering object (/steps_13_14_15_scenic_tf_regulon_analysis_scripts/L3_tfhlike_object_scenic) 
# 63 TFs from ATAC-based TF motif enrichment analysis of GNG4 RNA-defined groups - refer to code detailed above
# 2 TFs identified in common by each approach - NFATC1 and BACH1 - shown in schematic

# Next, the motifmatchr package was used to assess whether motifs for NFATC1 and BACH1 may be present specifically within the three GNG4 DAP regions of interest described above (PLS, pELS, and dELS-containing DAPs enriched in GC-like Tfh group)

# First, make GRanges object from positions of GC Tfh-enriched GNG4 DAPs
gng4_dap_strings <- c("chr1-235649534-235650662","chr1-235648709-235649522","chr1-235641919-235642869") # coordinates for PLS, pELS, and dELS-containing DAPs
gng4_dap_df <- do.call(rbind, strsplit(gng4_dap_strings, "-"))
colnames(gng4_dap_df) <- c("chr", "start", "end")
gng4_dap_df <- as.data.frame(gng4_dap_df)
gng4_dap_df$start <- as.numeric(gng4_dap_df$start)
gng4_dap_df$end <- as.numeric(gng4_dap_df$end)
gng4_dap_gr <- GRanges(seqnames = gng4_dap_df$chr, ranges = IRanges(start = gng4_dap_df$start, end = gng4_dap_df$end))

# Get TF motif PFM information from JASPAR2020
pfm_list <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
pfm_df <- do.call(rbind, lapply(pfm_list, function(m) { # creating dataframe of PFM IDs, common TF name, and TF family - 749 total from JASPAR2020
  data.frame(
    matrix_id = ID(m),
    tf_name = name(m),
    class = matrixClass(m),
    stringsAsFactors = FALSE
  )
}))

# Filter TF motif PFM information dataframe for NFATC1 and BACH1
gng4_tf <- c("NFATC1", "BACH1")
gng4_tf_ids <- pfm_df$matrix_id[pfm_df$tf_name %in% gng4_tf]
gng4_tf_motifs <- pfm_list[names(pfm_list) %in% gng4_tf_ids]
View(gng4_tf_motifs)
# NFATC1 - MA0624.1, Rel homology region (RHR) factors, NFAT-related factors
# BACH1 - MA1633.1, Basic leucine zipper factors (bZIP), Jun-related factors

# Score TF motifs for NFATC1 and BACH1 within GNG4 PLS, pELS, and dELS-containing DAP regions
gng4_dap_motif_matches <- matchMotifs(gng4_tf_motifs, gng4_dap_gr, genome = BSgenome.Hsapiens.UCSC.hg38, out = "scores", 
                                      p.cutoff = 0.05) # defining significance threshold as P < 0.05, permissive scan
rownames(gng4_dap_motif_matches) <- c('PLS','pELS','dELS')
colnames(gng4_dap_motif_matches) <- c("NFATC1", "BACH1")
gng4_dap_motif_matches@assays@data$motifMatches # matches found for all TF/DAP combinations, based on P < 0.05 cutoff
# Note P-values for all motif/DAP matches were < 0.001, as low as P < 0.0001 for the NFATC1/pELS match, and P < 0.00001 for the BACH1/PLS match
# Additional references for P-values used in motifmatchr analysis - https://github.com/GreenleafLab/motifmatchr/issues/3 and https://github.com/jhkorhonen/MOODS/issues/12#issuecomment-405912018
gng4_dap_motif_scores <- as.matrix(assay(gng4_dap_motif_matches)) # accessing motif scores assay by default, refer to gng4_dap_motif_matches@assays@data$motifScores

# Create color scheme for heatmap of motif scores in each DAP region
motif_max_val <- max(gng4_dap_motif_scores, na.rm = TRUE) # scaling color to maximum motif match score, from white/low to red/high
motif_heatmap_hi_col  <- brewer.pal(11, "RdBu")[4]
motif_heatmap_cols  <- colorRamp2(c(0, motif_max_val),
                                  c("white", motif_heatmap_hi_col))

# Assemble heatmap for motif match scores per TF and DAP region
pdf('/filepath/fig5/fig5c/gng4_tf_motifmatchr.pdf', width = 4.666, height = 4)
Heatmap(gng4_dap_motif_scores,
        col = motif_heatmap_cols,
        name = NULL,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(
          title = NULL,
          legend_height = unit(4, "cm"),
          labels_gp = gpar(fontsize = 12)
        ),
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.rect(x, y, width = w, height = h, # grid for heatmap
                    gp = gpar(col = "black", fill = NA, lwd = 1))
          grid.text(round(gng4_dap_motif_scores[i, j], 1), x, y, # motif scores from above
                    gp = gpar(fontsize = 20))
        })
dev.off()

# Fig 5D - L3 3WNN UMAP feature plots of SCENIC AUC scores and chromVAR deviation scores for NFATC1 and BACH1 ----

# Set color palette for feature plots
tf_umap_hi_col  <- brewer.pal(11, "RdBu")[2]
tf_umap_lo_col <- brewer.pal(11, "RdBu")[10]

# NFATC1 chromVAR UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'chromvar'
motif_scores <- FetchData(l3_teaseq_tfh_obj, vars = 'MA0624.1', slot = "data", assay = "chromvar")
min_cutoff <- quantile(motif_scores[[1]], probs = 0.01, na.rm = TRUE)
max_cutoff <- quantile(motif_scores[[1]], probs = 0.99, na.rm = TRUE)
pdf('/filepath/fig5/fig5d/nfatc1_chromvar_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'MA0624.1', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               min.cutoff = min_cutoff,
               max.cutoff = max_cutoff,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = tf_umap_lo_col, mid = "white", high = tf_umap_hi_col,
    midpoint = 0
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# BACH1 chromVAR UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'chromvar'
motif_scores <- FetchData(l3_teaseq_tfh_obj, vars = 'MA1633.1', slot = "data", assay = "chromvar")
min_cutoff <- quantile(motif_scores[[1]], probs = 0.01, na.rm = TRUE)
max_cutoff <- quantile(motif_scores[[1]], probs = 0.99, na.rm = TRUE)
pdf('/filepath/fig5/fig5d/bach1_chromvar_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'MA1633.1', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               min.cutoff = min_cutoff,
               max.cutoff = max_cutoff,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = tf_umap_lo_col, mid = "white", high = tf_umap_hi_col,
    midpoint = 0
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# Import SCENIC AUC score matrix, scale across cells for each TF regulon, and save as new assay for feature plots
DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
auc_matrix <- l3_teaseq_tfh_obj@assays$SCENIC@data # AUC scores range 0.0000000 to 0.8611048
scaled_auc <- t(scale(t(auc_matrix))) # range -3.238668 16.653320 (transpose matrix to get cell-by-regulon matrix, get regulon Z-score across cells, transpose back for format needed for Seurat)
l3_teaseq_tfh_obj[['scaled_scenic']] <- CreateAssayObject(data = scaled_auc)

# NFATC1 SCENIC UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'scaled_scenic'
pdf('/filepath/fig5/fig5d/nfatc1_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'NFATC1-extended (184g)', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = tf_umap_lo_col, mid = "white", high = tf_umap_hi_col,
    midpoint = 0
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# BACH1 SCENIC UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'scaled_scenic'
pdf('/filepath/fig5/fig5d/bach1_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'BACH1-extended (39g)', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = tf_umap_lo_col, mid = "white", high = tf_umap_hi_col,
    midpoint = 0
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# GNG4 RNA UMAP - scaled here, in contrast to in Fig 3E
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
pdf('/filepath/fig5/fig5d/gng4_rna_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'GNG4', order = TRUE, slot = 'scale.data'
,               legend.position = 'right', legend.length = 11,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = tf_umap_lo_col, mid = "white", high = tf_umap_hi_col,
    midpoint = 0
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# UMAP showing GC vs nonGC-like Tfh group annotations
pdf('/filepath/fig5/fig5d/gclike_class_dimplot.pdf', width = 6, height = 4)
Idents(l3_teaseq_tfh_obj) <- 'gc_vs_nongc_like'
do_DimPlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', group.by = 'gc_vs_nongc_like',
           colors.use = c('nonGC' = brewer.pal(11, "RdBu")[9],
                          'GC' = brewer.pal(10, "RdBu")[3]),
               border.size = 1.5, pt.size = 2.5) +
  coord_fixed() + NoAxes() + ggtitle(NULL) + NoLegend()
dev.off()

# Fig 5E - Experimental schematic for mononuclear cell in vitro stimulation ----

# PDF exported from Biorender

# Fig 5F - representative histograms of CD25, PD1, and Gy4 protein expression in gated CD4 T cell subsets with versus without stimulation ----

# PDF exported from FlowJo, further annotated in Illustrator

# Fig 5G - Barplot of % Gy4+ cells for each CD4 T cell subset in stimulated versus unstimulated conditions ----

# Import percentage values from FlowJo table export
tfh_gng4_flow_df <- read_excel("/filepath/fig5/fig5_flow_data_v2.xlsx", sheet = "Gy4_Stim_Freq")

# Format dataframe and subset names for barplot
tfh_gng4_flow_df_long <- tfh_gng4_flow_df %>%
  pivot_longer(cols = c(Naïve, nonTfh, Tfh),
               names_to = "Subset",
               values_to = "Frequency")

tfh_gng4_flow_df_long <- tfh_gng4_flow_df_long %>%
  mutate(
    Condition = factor(Condition, levels = c("Ctrl", "Stim")),
    Subset = factor(Subset, levels = c("Naïve", "nonTfh", "Tfh")),
    Tissue = factor(Tissue, levels = c("PBMC", "Tonsil")),
    Group = paste(Tissue, Subset, Condition, sep = "_")
  )

x_axis_order <- c(
  "PBMC_Naïve_Ctrl", "PBMC_Naïve_Stim",
  "PBMC_nonTfh_Ctrl", "PBMC_nonTfh_Stim",
  "PBMC_Tfh_Ctrl", "PBMC_Tfh_Stim",
  "Tonsil_Naïve_Ctrl", "Tonsil_Naïve_Stim",
  "Tonsil_nonTfh_Ctrl", "Tonsil_nonTfh_Stim",
  "Tonsil_Tfh_Ctrl", "Tonsil_Tfh_Stim"
)

tfh_gng4_flow_df_long <- tfh_gng4_flow_df_long %>%
  mutate(
    Group = paste(Tissue, Subset, Condition, sep = "_"),
    Group = factor(Group, levels = x_axis_order)
  )

# Assemble barplot
tfh_stim_gy4_freq_barplot <- ggplot(tfh_gng4_flow_df_long, aes(x = Group, y = Frequency, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar",
               position = position_dodge(width = 0.8),
               color = "black", width = 0.7) +
  geom_point(aes(group = Donor, fill = Condition),
             position = position_dodge(width = 0.8),
             shape = 21, color = "black", size = 3, stroke = 0.5) +
  geom_line(aes(group = interaction(Donor, Tissue, Subset)),
            position = position_dodge(width = 0.8),
            color = "black", linewidth = 0.4, alpha = 0.8) +
  scale_fill_manual(values = c("Ctrl" = "#7abdff", "Stim" = "#ff8686")) +
  scale_color_manual(values = c("Ctrl" = "#7abdff", "Stim" = "#ff8686")) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  theme_classic(base_size = 15, base_family = 'sans') +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + NoLegend()

tfh_stim_gy4_freq_barplot
pdf('/filepath/fig5/fig5g/tfh_stim_gy4_freq_barplot_2.pdf', width = 5.75, height = 4)
tfh_stim_gy4_freq_barplot
dev.off()

# Statistics for % Gy4+ differences between paired stimulated versus unstimulated conditions per donor were determined in Prism
# Paired t-test with two-stage step-up Benjamini, Krieger, and Yekutieli FDR correction for multiple comparisons, *Q < 0.05, **Q < 0.01
# *Q and **Q symbols added in Illustrator

# Gy4 GMFI values for each subset were analyzed and visualized in Prism, shown in Fig. S22E (tonsil) and S22F (PBMC)

# Fig 5H - GNG4 ATAC coverage plot for L4 GC vs nonGC-like Tfh groups, with rheumatoid arthritis fine-mapping GWAS variant positions annotated, GNG4 DAP regions highlighted, and co-accessible regions indicated by Cicero links ----

# Find networks of co-accessible peaks using Cicero on L4 Tfh-only object without Tcm
l4_teaseq_tfh_no_tcm_obj.cds <- as.cell_data_set(l4_teaseq_tfh_no_tcm_obj) # convert to cds object to run Cicero
tfh_no_tcm_cicero <- make_cicero_cds(l4_teaseq_tfh_no_tcm_obj.cds, reduced_coordinates = reducedDims(l4_teaseq_tfh_no_tcm_obj.cds)$UMAP.WNN.HARMONY) # use 3WNN UMAP embedding
sl <- seqlengths(Hsapiens)
genome.df <- data.frame(
  chr    = names(sl),
  length = as.numeric(sl),
  stringsAsFactors = FALSE
)
tfh_no_tcm_cicero_conns <- run_cicero(tfh_no_tcm_cicero, genomic_coords = genome.df, sample_num = 100)
tfh_no_tcm_cicero_ccans <- generate_ccans(tfh_no_tcm_cicero_conns)
tfh_no_tcm_cicero_links <- ConnectionsToLinks(conns = tfh_no_tcm_cicero_conns, ccans = tfh_no_tcm_cicero_ccans)
Links(l4_teaseq_tfh_no_tcm_obj) <- tfh_no_tcm_cicero_links # replace 'Links' slot with CCAN links computed from Cicero analysis

# Find Cicero CCAN links falling within the GNG4 locus of L4 Tfh-only cells
ann <- Annotation(l4_teaseq_tfh_no_tcm_obj)
gng4_range <- ann[ann$gene_name == "GNG4"] # get GNG4 locus coordinates
all_links <- Links(l4_teaseq_tfh_no_tcm_obj) #  get coordinates for all links
gng4_links <- subsetByOverlaps(all_links, gng4_range, ignore.strand = TRUE) # get links falling within the GNG4 locus
gng4_links_df <- as.data.frame(gng4_links)

# Save Cicero CCAN links for Data File S7
gng4_links_df$strand <- NULL
gng4_links_df_save <- gng4_links_df %>%
  rename(
    chromosome = seqnames,
    cicero_peak_pair_coaccess_score = score,
    link_start = start,
    link_end = end,
    link_width = width,
    cis_co_access_network_group_id = group
  )
gng4_links_df_save <- gng4_links_df_save %>% arrange(desc(cicero_peak_pair_coaccess_score))
View(gng4_links_df_save)
saveRDS(gng4_links_df_save, "/filepath/fig5/fig5h_l4_tfh_cicero_ccans_gng4_link_stats.rds")

# Compile Cicero CCAN information and all other GNG4 gene regulation files above into tabs of single spreadsheet for Data File S7
df_dir   <- '/filepath/fig5/'
df_files <- list.files(df_dir, pattern = "\\.rds$", full.names = TRUE)
xlsx_sheets <- setNames(vector("list", length(df_files)), nm = basename(df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(df_files)) {
  df <- readRDS(df_files[i])
  xlsx_sheets[[i]] <- df
}
write_xlsx(xlsx_sheets, path = "/filepath/fig5/fig5_tfh_gng4_gene_reg_dap_links_scenic_cicero_spreadsheets.xlsx")  

# Finding all ATAC peaks within GNG4 region
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ATAC'
tfh_annotation <- Annotation(l4_teaseq_tfh_no_tcm_obj)
gng4_range <- tfh_annotation[tfh_annotation$gene_name == 'GNG4', ]
tfh_peak_ranges <- granges(l4_teaseq_tfh_no_tcm_obj[["ATAC"]]) # get GRanges for all peaks within ATAC assay
tfh_peak_gng4 <- subsetByOverlaps(tfh_peak_ranges, gng4_range) # filter for peaks within GNG4
tfh_peak_gng4_df <- as.data.frame(tfh_peak_gng4)
View(tfh_peak_gng4_df) # 14 total peak regions in GNG4 identified by cellranger-arc
gc_vs_nongc_tfh_dap_save %>% filter(p_val_adj < 0.05 & closest_gene == 'GNG4') %>% View() # of these 14, only 3 were found to be significantly enriched in the GC-like L4 Tfh group (refer to Fig 5A analysis above)

# Additional peaks of interest further downstream of the GNG4 promoter region that were not noted in Fig 5A 
# chr1-235567556-235568514, far beyond region with significant accessibility in L4 Tfh, not included in visualization window
# chr1-235606974-235608281, far beyond region with significant accessibility in L4 Tfh, not included in visualization window
# chr1-235632708-235633527, dELS-containing but low accessibility in L4 Tfh
# chr1-235635874-235637036, dELS-containing but low accessibility in L4 Tfh

# Create GRanges object for GNG4 DAPs of interest for visualization
gng4_dap_df <- data.frame(
  seqnames = "chr1",
  start = c(235641919, 235649534, 235648709, 235632708, 235635874), # start coordinates for dELS, PLS, pELS, low signal peak #1, low signal peak #2
  end   = c(235642869, 235650662, 235649522, 235633527, 235637036) # end coordinates for dELS, PLS, pELS, low signal peak #1, low signal peak #2
)
gng4_dap_gr <- makeGRangesFromDataFrame(gng4_dap_df)
rdbu_colors <- brewer.pal(11, "RdBu")
gng4_dap_gr$color <- c('darkgrey',rdbu_colors[3],'darkgrey','darkgrey','darkgrey') # coloring PLS red, other DAP regions grey for visualization

# Having identified co-accessible and differentially accessible regions of GNG4 in L4 GC-like Tfh, next we considered whether any known variants from GWAS/eQTL studies are found near these regions of interest

# Ishigaki, Sakaue, Terao et al. Multi-ancestry genome-wide association analyses identify novel genetic mechanisms in rheumatoid arthritis. Nat. Genet. 2022. PMID 36333501. https://pubmed.ncbi.nlm.nih.gov/36333501/
# Fine-mapping GWAS results include set of variants mapped to GNG4 that are associated with RA diagnosis
# Lead variant annotated as rs1188620266 (chr1:235800357:CAA:C) in original study
# Due to changes in gnomAD versions since this study and our reanalysis, we have chosen to refer to this variant using the updated rs61512163 identifier (chr1:235637057:CAA:C), which we confirmed with the original study authors
# Refer to Supplemental Materials and Data File S9 for additional GWAS and variant details

# rsIDs for GNG4 variants identified in RA GWAS, where available
gng4_snps <- c(
  'rs61512163','rs12133886','rs35458456','rs34999365', # note rs61512163 rsID used here rather than rs1188620266
  'rs10926320','N/A','rs10802904','rs6429213',
  'rs12133526','N/A','N/A','rs7555242',
  'rs4391655'
)

# Positions of GNG4 variants in hg38, determined using Broad Institute Liftover Webtool (refer to Data File S9)
gng4_snp_pos <- c(
  235637057, 235638484, 235638486, 235638482,
  235638201, 235637653, 235643279, 235642088,
  235643355, 235635572, 235631575, 235635266,
  235637611
)

# rsID/positions were visualized in L4 Tfh GNG4 locus coverage plot, then further annotated as below
# Within GNG4 dELS-containing DAP region - rs6429213 (235642088)
# Proximal to same GNG4 dELS-containing DAP region - rs10802904 (235643279) and rs12133526 (235643355)

# Assigning colors to each variant for exploration - each thin colored line was then masked in Illustrator where we used thicker gold lines
gng4_snp_gr_colors <- c('green4','darkorange3','darkorange3','darkorange3',
                        'darkorange3','purple3','red','blue',
                        'cyan','purple3','purple3','darkorange3',
                        'darkorange3')
# green4 = lead variant reported in study
# purple3 = variants lacking rsIDs
# red = rs10802904, proximal to GNG4 dELS-containing DAP region
# blue = rs6429213, within dELS-containing DAP region
# cyan = rs12133526, proximal to GNG4 dELS-containing DAP region
# darkorange3 = all other variants

# Make GRanges object for RA GWAS variants
gng4_snp_gr <- GRanges(
  seqnames = Rle("chr1", length(gng4_snps)),
  ranges   = IRanges(start = gng4_snp_pos - 10, end = gng4_snp_pos + 10), # widening single nucleotide positions to make colored lines visible
  snp      = gng4_snps
)
gng4_snp_gr$color <- gng4_snp_gr_colors

# Combine GRanges objects for RA GWAS variants and GNG4 L4 GC-like Tfh DAPs
combined_gr <- c(gng4_snp_gr,gng4_dap_gr)

# Determine visualization window for coverage plot
start_region <- min(start(combined_gr)) - 2000
end_region   <- max(start(combined_gr)) + 2000
plot_region <- paste0("chr1-", start_region, "-", end_region)

# Assemble coverage plot
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'ATAC'
l4_teaseq_tfh_no_tcm_obj$gc_vs_nongc_like <- factor(l4_teaseq_tfh_no_tcm_obj$gc_vs_nongc_like,  levels = c('GC','nonGC'))
gng4_snp_plot <- CoveragePlot(l4_teaseq_tfh_no_tcm_obj, region = plot_region, links = TRUE, annotation = TRUE, peaks = FALSE, layer = 'data', features = NULL,
                              group.by = 'gc_vs_nongc_like', region.highlight = combined_gr) & scale_fill_manual(values = c('GC' = 'black', 'nonGC' = 'black')) 
gng4_snp_plot[[1]][[1]] <- gng4_snp_plot[[1]][[1]] + theme(axis.ticks = element_blank()) + labs(y = 'Normalized Peak Signal (0-230)')
gng4_snp_plot[[1]][[2]] <- gng4_snp_plot[[1]][[2]] + scale_color_manual(values = c('black','black'))
gng4_snp_plot_links <- gng4_snp_plot[[1]][[3]]
gng4_snp_plot_links <- gng4_snp_plot_links +
  scale_color_gradient(
    low  = "lightgrey",
    high = rdbu_colors[3],
    name = "Cicero Score" 
  ) +
  theme(
    legend.key.size   = unit(0.3, "cm"),
    axis.text.x = element_text(size = 5),
    axis.title.x = element_blank(),
    legend.text       = element_text(size = 8),
    legend.title      = element_blank()
  )
gng4_snp_plot[[1]][[3]] <- gng4_snp_plot_links
pdf('/filepath/fig5/fig5h/dao_gng4_eqtl_ra_gwas_snp_ld_reanalysis.pdf', width = 5, height = 4)
gng4_snp_plot
dev.off()

# Data File S9 shows linkage disequilibrium analysis between these RA GWAS variants and the eQTL variants found in reanalysis of data from Schmiedel et al. Sci Immunol 2022.
# An additional variant of GNG4 (rs551907784) was reported in Verma et al. Science 2024 GWAS, but not found to be in LD with the eQTL variants, and was not visualized

# Fig 5I - Reanalysis of longitudinal SARS2 mRNA vaccination LN FNA scRNAseq study ----

# Reanalysis of Borcherding et al. CD4+ T cells exhibit distinct transcriptional phenotypes in the lymph nodes and blood following mRNA vaccination in humans. Nat Immunol. 2024. PMID 39164479. https://pubmed.ncbi.nlm.nih.gov/39164479/
# Import T cell subclustering Seurat object from SARS2 mRNA vaccination LN FNA scRNAseq study
# .rds file obtained from Figure2 tab of https://cellpilot.emed.wustl.edu/
covid_tcell_obj <- readRDS('/filepath/reanalysis/borcherding_mrna_vax/BorcherdingFig2.rds')

# Refer to Fig S24 and related code in 'figS24_reanalysis_human_vaccination_scrnaseq_gng4_gctfh.R' for further annotation details of GC Tfh cluster

# Author annotations list c3 as CD4+ GC Tfh - our reanalysis supports this annotation, as shown by dotplot below
Idents(covid_tcell_obj) <- 'seurat_clusters' # pre-existing clusters c0-c11
DotPlot(covid_tcell_obj, features = c('GNG4','TIGIT','PDCD1','CXCR5','TOX2','IL21','ASCL2','B3GAT1','NFATC1','BACH1')) # GC signature also observed to lesser extent in c11_TfhPro and c8_IL10Tfh

# Here we group all nonGC Tfh clusters into 'Other' for comparison with the GCTfh cluster over time.
covid_tcell_obj <- SetIdent(
  covid_tcell_obj,
  value = ifelse(Idents(covid_tcell_obj) == '3', "GCTfh", "Other")
)
covid_tcell_obj$gctfh_vs_other <- Idents(covid_tcell_obj)
table(covid_tcell_obj$gctfh_vs_other)

# In our reanalysis, we found that the d110 timepoint featured an outlier low number of GCTfh, potentially reflecting poor LN sampling during the FNA or other timepoint-specific sample issues.
# Moreover, the authors excluded the d110 timepoint in Fig 4 analysis 'dLN CD4+ TFH cells from day 110 after vaccination included significantly fewer spike-specific cells (10) than cells from days 28 (94), 60 (64) and 201 (70) and were therefore excluded from the analysis.'
covid_gctfh_freq_df <- covid_tcell_obj@meta.data %>%
  group_by(timepoint) %>%
  summarise(
    n_GCTfh = sum(gctfh_vs_other == "GCTfh"),   
    n_Other = sum(gctfh_vs_other == "Other"), 
    percent_GCTfh = n_GCTfh / (n_GCTfh + n_Other) * 100
  ) %>%
  mutate(
    cell_tot = n_GCTfh + n_Other
  )
covid_gctfh_freq_df # d110 has 3.33% GCTfh
IQR(covid_gctfh_freq_df$percent_GCTfh) # % GCTfh IQR = 3.93859
quantile(covid_gctfh_freq_df$percent_GCTfh, 0.25) - 1.5*IQR(covid_gctfh_freq_df$percent_GCTfh) # 4.274979 % GCTFh lower bound outlier threshold. 3.33% at d110 < 4.274979% lower bound threshold, suggesting that d110 may be an outlier

# Append timepoint to GCTfh/Other annotation
covid_tcell_obj$cluster_timepoint <- paste0(
  covid_tcell_obj$gctfh_vs_other,
  "_",
  covid_tcell_obj$timepoint
)

# Given our outlier analysis and analyses performed in the original study, we filtered out d110 from our downstream reanalysis 
covid_tcell_obj_filt <-  subset(covid_tcell_obj, timepoint != 'd110')

# Create dotplot of RNA expression for GNG4 and other features in GCTfh vs nonGCTfh and nonTfh over time
covid_tcell_obj_filt$cluster_timepoint <- factor(covid_tcell_obj_filt$cluster_timepoint, levels = rev(c('GCTfh_d201','GCTfh_d60','GCTfh_d28',
                                                                                                        'Other_d201','Other_d60','Other_d28')))
covid_tcell_labels <- sub(".*_d", "", levels(covid_tcell_obj_filt$cluster_timepoint))
rdbu_colors <- brewer.pal(11, "RdBu")
borcherding_dotplot_feats <- c('IL21','TOX2','BCL6','S1PR2','B3GAT1','NFATC1','BACH1','GNG4') # features of interest
pdf('/filepath/fig5/fig5i/borcherding_covid_fna_reanalysis_dotplot.pdf', width = 6.25, height = 6.5)
borcherding_covid_fna_reanalysis_dotplot <- DotPlot(covid_tcell_obj_filt, features = borcherding_dotplot_feats, group.by = 'cluster_timepoint', scale.by = 'size', dot.scale = 12, cols = 'RdBu') +
  theme(legend.justification.right = 'bottom',
        legend.title = element_text(size = 12, family = 'sans'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 14, face = 'italic', family = 'sans', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, family = 'sans')) +
  scale_y_discrete(labels = covid_tcell_labels) +
  scale_color_gradient2(
    low = rdbu_colors[11],           
    mid = "white",                
    high = rdbu_colors[2],        
    midpoint = 0                  
  )
borcherding_covid_fna_reanalysis_dotplot + guides( # adjusting legend
  size = guide_legend(
    title = "% Expr",
    order = 1,
    title.position = "top",
    title.hjust    = 0.5,
    reverse = TRUE
  ),
  color = guide_colorbar(
    title = "Scaled Expr",
    order = 2,
    title.position = "top",
    title.hjust    = 0.5
  )
) 
dev.off()

# Determining fold change in % Expr and Avg Expr for GNG4 in GCTfh from d28 to d60 and from d60 to d201
# Get percent and average expression values output by DotPlot function
borcherding_gctfh_time_bin_stats <- DotPlot(covid_tcell_obj_filt,
                                            features = borcherding_dotplot_feats, 
                                            group.by = 'cluster_timepoint',
                                            scale = FALSE) # using normalized expression values (avg.exp), not Z-scores column (avg.exp.scaled)
borcherding_gctfh_time_bin_stats$data

# d60 / d28 Pct Expr
36.3225806/27.2922636 # 1.330875 fold increase
# d201 / d60 Pct Expr
43.8006952/36.3225806 # 1.205881 fold increase

# d60 / d28 Avg Expr
1.59654083/1.28730403 # 1.24022 fold increase
# d201 / d60 Avg Expr
1.96361706/1.59654083 # 1.22992 fold increase

# Fig 5J - Reanalysis of quadrivalent inactivated influenza longitudinal FNA scRNAseq data ----

# Reanalysis of Schattgen et al. Influenza vaccination stimulates maturation of the human T follicular helper cell response. Nat Immunol. 2024. PMID 39164477. https://pubmed.ncbi.nlm.nih.gov/39164477/
# Import Seurat objects for T cell clustering and Level 2 Tfh subclustering from IIV LN FNA scRNAseq study 
# Seurat object .rds files ('intergrated_Tcells_harmony_bydonor.rds' and 'intergrated_Tfh_harmony_bydonor.rds') obtained from Zenodo data repository (Version 3, Oct 20 2023, https://zenodo.org/records/12611325)
flu_tcell_obj <- readRDS('/filepath/reanalysis/intergrated_Tcells_harmony_bydonor.rds')
flu_tfh_obj <- readRDS('/filepath/reanalysis/intergrated_Tfh_harmony_bydonor.rds')

# Note - as provided, the Level 1 T cell object does not have a specifically annotated 'GC Tfh' cluster, whereas the Level 2 Tfh subclustering object does ('GC')
# Therefore, in steps below we examine expression of GNG4 and other GC-like Tfh features in both the L2 (S24D-F) and L1 objects (S24G-I), then map the 'GC' cluster from L2 to L1 (S24J-K)
# Refer to Fig S24 and related code in 'figS24_reanalysis_human_vaccination_scrnaseq_gng4_gctfh.R' for further annotation details of GC Tfh cluster

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

# Get number and proportion of L2-to-L1 mapped 'GC Tfh' vs other clusters in L1 object over time
flu_tfh_freq_df <- flu_tcell_obj@meta.data %>%
  group_by(time) %>%
  summarise(
    n_GCTfh = sum(gctfh_annot == "yes"),   
    n_Other_Cells = sum(gctfh_annot == "no"), 
    percent_GCTfh = n_GCTfh / (n_GCTfh + n_Other_Cells) * 100
  )
flu_tfh_freq_df <- flu_tfh_freq_df %>%
  mutate(cell_tot = n_GCTfh + n_Other_Cells)
flu_tfh_freq_df <- flu_tfh_freq_df %>% arrange(cell_tot)
flu_tfh_freq_df # y1_d26 features only 722 cells, of which only 7 were GC Tfh - excluding this timepoint from downstream analyses of RNA expression over time due to low sampling

# Remove low cell number d26 timepoint from object
Idents(flu_tcell_obj) <- 'day'
flu_tcell_obj_filt <- subset(flu_tcell_obj, day != '26')

# Create GCTfh or nonGCTfh_nonTfh vs timepoint identifier
flu_tcell_obj_filt$gctfh_yesno_day <- paste0(
  flu_tcell_obj_filt$gctfh_annot, # GC Tfh cell identity determined above
  "_",
  flu_tcell_obj_filt$day # merging sampling days if shared between the two years of the study
)
table(flu_tcell_obj_filt$gctfh_yesno_day, useNA = "ifany")

# Set GC Tfh yes/no_day# identity order for dotplot
flu_tcell_obj_filt$gctfh_yesno_day <- factor(flu_tcell_obj_filt$gctfh_yesno_day, levels = 
                                                  rev(c('yes_180','yes_120','yes_90','yes_60','yes_28','yes_14','yes_12','yes_7','yes_5','yes_0',
                                                        'no_180','no_120','no_90','no_60','no_28','no_14','no_12','no_7','no_5','no_0')))
levels(flu_tcell_obj_filt$gctfh_yesno_day)

# Make abbreviated labels for dotplot
short_labels <- sub(".*_", "", levels(flu_tcell_obj_filt$gctfh_yesno_day))

# Set color palette for dotplot
rdbu_colors <- brewer.pal(11, "RdBu")

# Select RNA features to visualize
schattgen_dotplot_feats <- c('IL21','TOX2','BCL6','S1PR2','B3GAT1','NFATC1','BACH1','GNG4')

# Create dotplot
flu_fna_reanalysis_dotplot <- DotPlot(flu_tcell_obj_filt, features = schattgen_dotplot_feats, group.by = 'gctfh_yesno_day', scale.by = 'size', dot.scale = 9, cols = 'RdBu') +
  theme(legend.justification.right = 'bottom',
        legend.title = element_text(size = 12, family = 'sans'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 14, face = 'italic', family = 'sans', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, family = 'sans')) +
  scale_y_discrete(labels = short_labels) +
  scale_color_gradient2(
    low = rdbu_colors[11],
    mid = "white",
    high = rdbu_colors[2],
    midpoint = 0
  )

# Adjust legends for dotplot - further annotated in Illustrator
pdf('/filepath/fig5/fig5j/flu_fna_reanalysis_dotplot_smallerdots.pdf', width = 6.25, height = 6.5)
flu_fna_reanalysis_dotplot + guides(
  size = guide_legend(
    title = "% Expr",
    order = 1,
    title.position = "top",
    title.hjust    = 0.5,
    reverse = TRUE
  ),
  color = guide_colorbar(
    title = "Scaled Expr",
    order = 2,
    title.position = "top",
    title.hjust    = 0.5
  )
) 
dev.off()

# Calculate differences in Avg and Pct Expr of GNG4 at d0-28 vs d60+ timepoints

# Bin d0-28 vs d60+ timepoints as 'Early' vs 'Late' in GC response for comparison
flu_tcell_obj_filt$time_bin <- ifelse(
  as.numeric(sub(".*_", "", flu_tcell_obj_filt$gctfh_yesno_day)) <= 28,
  "0-28", "60+"
)
table(flu_tcell_obj_filt$time_bin)

# Subset on GC Tfh to calculate differences in GNG4 RNA expression over time
flu_tcell_obj_filt_tfh_subset <- subset(flu_tcell_obj_filt, gctfh_annot == 'yes')

# Factor time bins for comparison
flu_tcell_obj_filt_tfh_subset$time_bin <- factor(
  flu_tcell_obj_filt_tfh_subset$time_bin,
  levels = c("0-28","60+")
)

# Get percent and average expression values output by DotPlot function
schattgen_gctfh_time_bin_stats <- DotPlot(
  object  = flu_tcell_obj_filt_tfh_subset,
  features = schattgen_dotplot_feats,
  group.by = "time_bin",
  scale = FALSE # using normalized expression values (avg.exp), not Z-scores column (avg.exp.scaled)
)

# Determine differences and add in Illustrator to dotplot
schattgen_gctfh_time_bin_stats$data
# GNG4 Avg Expr d0-28 vs d60+ = 0.26100881 vs 1.17268623 (4.4928990328x)
# GNG4 % Expr d0-28 vs d60+ = 11.769416 vs 42.776524 (3.6345494118x)