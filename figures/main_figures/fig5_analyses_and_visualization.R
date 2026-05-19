# Sam Barnett Dubensky et al.
# Derek A. Oldridge & Laura A. Vella Labs at the Children's Hospital of Philadelphia
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned Tfh in human lymphoid tissue
# Code and data visualization for Fig. 5
# Fig. 5 - GNG4 expression delineates a highly activated state within the GC-like Tfh pool

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
library(writexl) # 1.5.4
library(ggnewscale) # 0.5.0
library(readxl) # 1.4.3
library(UCell) # 2.8.0
library(scales) # 1.3.0
library(rlang) # 1.1.4
library(ggrepel) # 0.9.5
library(speckle) # 1.4.0

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
# Refer to Methods text for DAP region annotation - visualized DAP regions in hg38 using UCSC Genome Browser with track for ENCODE Candidate Cis-Regulatory Elements (cCREs across human cell types) 
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

# UMAP showing GNG4 RNA+ vs RNA- states
DefaultAssay(l3_teaseq_tfh_obj) <- 'RNA'
l3_teaseq_tfh_obj$gng4_rna_class <- ifelse(FetchData(l3_teaseq_tfh_obj, vars = "rna_GNG4") >= 1, "GNG4_RNA_high", "GNG4_RNA_low")
Idents(l3_teaseq_tfh_obj) <- 'gng4_rna_class'
pdf('/filepath/fig5/fig5c/gng4rna_class_dimplot.pdf', width = 6, height = 4)
do_DimPlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', group.by = 'gng4_rna_class', order = c('GNG4_RNA_high'), shuffle = FALSE,
           colors.use = c('GNG4_RNA_low' = 'lightgrey', 'GNG4_RNA_high' = '#d5604d'),
           border.size = 1.5, pt.size = 2.5) +
  coord_fixed() + NoAxes() + ggtitle(NULL) + NoLegend()
dev.off()

# Color scheme for chromVAR and SCENIC UMAP visualization
umap_hi_col  <- brewer.pal(11, "RdBu")[2]
umap_lo_col <- brewer.pal(11, "RdBu")[10]

# NFATC1 chromVAR UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'chromvar'
motif_scores <- FetchData(l3_teaseq_tfh_obj, vars = 'MA0624.1', slot = "data", assay = "chromvar")
min_cutoff <- quantile(motif_scores[[1]], probs = 0.01, na.rm = TRUE)
max_cutoff <- quantile(motif_scores[[1]], probs = 0.99, na.rm = TRUE)
pdf('/filepath/fig5/fig5c/nfatc1_chromvar_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'MA0624.1', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               min.cutoff = min_cutoff,
               max.cutoff = max_cutoff,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = umap_lo_col, mid = "white", high = umap_hi_col,
    midpoint = 0,
    labels = scales::label_number(accuracy = 0.1)
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
pdf('/filepath/fig5/fig5c/bach1_chromvar_umap.pdf', width = 6, height = 4)
motif_scores <- FetchData(l3_teaseq_tfh_obj, vars = 'MA1633.1', slot = "data", assay = "chromvar")
min_cutoff <- quantile(motif_scores[[1]], probs = 0.01, na.rm = TRUE)
max_cutoff <- quantile(motif_scores[[1]], probs = 0.99, na.rm = TRUE)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'MA1633.1', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               min.cutoff = min_cutoff,
               max.cutoff = max_cutoff,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = umap_lo_col, mid = "white", high = umap_hi_col,
    midpoint = 0,
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# NFATC1 SCENIC UMAP
DefaultAssay(l3_teaseq_tfh_obj) <- 'SCENIC'
auc_matrix <- l3_teaseq_tfh_obj@assays$SCENIC@data 
scaled_auc <- t(scale(t(auc_matrix)))
l3_teaseq_tfh_obj[['scaled_scenic']] <- CreateAssayObject(data = scaled_auc)
DefaultAssay(l3_teaseq_tfh_obj) <- 'scaled_scenic'
pdf('/filepath/fig5/fig5c/nfatc1_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'NFATC1-extended (184g)', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = umap_lo_col, mid = "white", high = umap_hi_col,
    midpoint = 0,
    labels = scales::label_number(accuracy = 0.1)
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
pdf('/filepath/fig5/fig5c/bach1_scenic_umap.pdf', width = 6, height = 4)
do_FeaturePlot(l3_teaseq_tfh_obj, reduction = 'umap.wnn.harmony', features = 'BACH1-extended (39g)', order = TRUE, 
               legend.position = 'right', legend.length = 11,
               border.size = 1.5, pt.size = 2.5) +
  scale_color_gradient2(
    low = umap_lo_col, mid = "white", high = umap_hi_col,
    midpoint = 0,
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme(
    legend.position = c(1, 0.5),
    legend.direction = 'vertical',
    plot.margin = margin(r = 10)
  ) +
  coord_fixed() + NoAxes() + ggtitle(NULL)
dev.off()

# Compile GNG4 gene regulation files above into tabs of single spreadsheet for Data File S7
df_dir   <- '/filepath/fig5/'
df_files <- list.files(df_dir, pattern = "\\.rds$", full.names = TRUE)
xlsx_sheets <- setNames(vector("list", length(df_files)), nm = basename(df_files) %>% str_remove("\\.rds$"))
for (i in seq_along(df_files)) {
  df <- readRDS(df_files[i])
  xlsx_sheets[[i]] <- df
}
write_xlsx(xlsx_sheets, path = "/filepath/fig5/fig5_tfh_gng4_gene_reg_dap_links_scenic_spreadsheets.xlsx")  

# Fig 5D - paired donor-level comparison of NFATC1 and BACH1 chromVAR scores in GNG4- versus GNG4+ subsets of tonsillar GC-like Tfh ----

# Extract tonsillar GC-like Tfh from L4 3WNN Tfh-only Seurat object
tonsil_gctfh_obj <- subset(l4_teaseq_tfh_no_tcm_obj, subset = hto.tissue == 'Tonsil' & gc_vs_nongc_like == 'GC')

# Specify colors for GNG4 RNA+ vs RNA- cells
gng4_class_cols <- c(
  "GNG4_RNA_low" = "lightgrey",
  "GNG4_RNA_high" = "#d5604d"
)

# Map JASPAR PFM IDs in chromVAR count matrix to motif names for BACH1 and NFATC1
motif_id_map <- c(
  "MA1633.1" = "BACH1",
  "MA0624.1" = "NFATC1"
)

# Get chromVAR scores for NFATC1/BACH1 and cell metadata 
motifs_use <- names(motif_id_map)
DefaultAssay(tonsil_gctfh_obj) <- "chromvar"
chromvar_df <- FetchData(
  tonsil_gctfh_obj,
  vars = c(motifs_use, "hto.donor", "gng4_rna_class"),
  slot = "data"
) %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    gng4_rna_class = factor(
      gng4_rna_class,
      levels = c("GNG4_RNA_low", "GNG4_RNA_high")
    )
  )

# Get mean chromVAR motif scores for each tonsil donor and GNG4 class
donor_chromvar_df <- chromvar_df %>%
  tidyr::pivot_longer(
    cols = all_of(motifs_use),
    names_to = "motif_id",
    values_to = "chromvar_deviation"
  ) %>%
  dplyr::mutate(
    motif = dplyr::recode(motif_id, !!!motif_id_map),
    motif = factor(motif, levels = c("NFATC1", "BACH1"))
  ) %>%
  dplyr::group_by(
    hto.donor,
    gng4_rna_class,
    motif
  ) %>%
  dplyr::summarise(
    mean_chromvar = mean(chromvar_deviation),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

# Assemble scatter plot
pdf('/filepath/fig5/fig5d/donor_chromvar_motif_scores_wide.pdf',
    width = 8,
    height = 3)
ggplot(
  donor_chromvar_df,
  aes(x = gng4_rna_class, y = mean_chromvar)
) +
  geom_line(
    aes(group = hto.donor),
    color = "grey",
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_point(
    aes(fill = gng4_rna_class),
    shape = 21,
    color = "black",
    stroke = 0.7,
    size = 8
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.3,
    linewidth = 0.3,
    alpha = 0.8,
    color = "grey"
  ) +
  facet_wrap(~ motif, scales = "free_y") +
  scale_fill_manual(values = gng4_class_cols) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  theme_classic() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    text = element_text(size = 18, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color = "black"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
dev.off()

# Prepare table to perform paired donor-level t-test
gng4_wide_df <- donor_chromvar_df %>%
  dplyr::select(
    hto.donor,
    motif,
    gng4_rna_class,
    mean_chromvar
  ) %>%
  tidyr::pivot_wider(
    names_from = gng4_rna_class,
    values_from = mean_chromvar
  )

# Perform paired donor-level t-test with BH FDR correction
chromvar_tests_gng4 <- gng4_wide_df %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(
    contrast = "GNG4_RNA_high vs GNG4_RNA_low",
    n_pairs = dplyr::n(),
    mean_GNG4_RNA_low = mean(GNG4_RNA_low),
    mean_GNG4_RNA_high = mean(GNG4_RNA_high),
    mean_paired_difference = mean(GNG4_RNA_high - GNG4_RNA_low),
    sd_paired_difference = sd(GNG4_RNA_high - GNG4_RNA_low),
    p_value = ifelse(
      n_pairs >= 3,
      t.test(
        x = GNG4_RNA_high,
        y = GNG4_RNA_low,
        paired = TRUE
      )$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_adj_BH)

# Save results
write_xlsx(donor_chromvar_df, path = '/filepath/fig5/fig5d/tonsil_gcliketfh_nfatc1_bach1_chromvar_donor_scores.xlsx')
write_xlsx(chromvar_tests_gng4, path = '/filepath/fig5/fig5d/tonsil_gcliketfh_nfatc1_bach1_chromvar_donor_statistics.xlsx')

# Fig 5D - paired donor-level comparison of NFATC1 and BACH1 SCENIC regulon AUC scores in GNG4- versus GNG4+ subsets of tonsillar GC-like Tfh ----

# Map regulon names in cell-by-SCENIC AUC score matrix to motif names for BACH1 and NFATC1
regulon_id_map <- c(
  "BACH1-extended (39g)" = "BACH1",
  "NFATC1-extended (184g)" = "NFATC1"
)

# Get SCENIC scores for NFATC1/BACH1 and cell metadata 
motifs_use <- names(regulon_id_map)
DefaultAssay(tonsil_gctfh_obj) <- "SCENIC"
scenic_df <- FetchData(
  tonsil_gctfh_obj,
  vars = c(motifs_use, "hto.donor", "gng4_rna_class"),
  slot = "data"
) %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    gng4_rna_class = factor(
      gng4_rna_class,
      levels = c("GNG4_RNA_low", "GNG4_RNA_high")
    )
  )

# Get mean SCENIC regulon AUC scores for each tonsil donor and GNG4 class
donor_scenic_df <- scenic_df %>%
  tidyr::pivot_longer(
    cols = all_of(motifs_use),
    names_to = "motif_id",
    values_to = "scenic_auc"
  ) %>%
  dplyr::mutate(
    motif = dplyr::recode(motif_id, !!!regulon_id_map),
    motif = factor(motif, levels = c("NFATC1", "BACH1"))
  ) %>%
  dplyr::group_by(
    hto.donor,
    gng4_rna_class,
    motif
  ) %>%
  dplyr::summarise(
    mean_scenic_auc = mean(scenic_auc),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

# Assemble scatter plot
pdf('/filepath/fig5/fig5d/donor_scenic_auc_scores_wide.pdf',
    width = 8,
    height = 3)

ggplot(
  donor_scenic_df,
  aes(x = gng4_rna_class, y = mean_scenic_auc)
) +
  geom_line(
    aes(group = hto.donor),
    color = "grey",
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_point(
    aes(fill = gng4_rna_class),
    shape = 21,
    color = "black",
    stroke = 0.7,
    size = 8
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.3,
    linewidth = 0.3,
    alpha = 0.8,
    color = "grey"
  ) +
  facet_wrap(~ motif, scales = "free_y") +
  scale_fill_manual(values = gng4_class_cols) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme_classic() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    text = element_text(size = 18, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color = "black"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )
dev.off()

# Prepare table to perform paired donor-level t-test
gng4_wide_df <- donor_scenic_df %>%
  dplyr::select(
    hto.donor,
    motif,
    gng4_rna_class,
    mean_scenic_auc
  ) %>%
  tidyr::pivot_wider(
    names_from = gng4_rna_class,
    values_from = mean_scenic_auc
  )

# Perform paired donor-level t-test with BH FDR correction
scenic_tests_gng4 <- gng4_wide_df %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(
    contrast = "GNG4_RNA_high vs GNG4_RNA_low",
    n_pairs = dplyr::n(),
    mean_GNG4_RNA_low = mean(GNG4_RNA_low),
    mean_GNG4_RNA_high = mean(GNG4_RNA_high),
    mean_paired_difference = mean(GNG4_RNA_high - GNG4_RNA_low),
    sd_paired_difference = sd(GNG4_RNA_high - GNG4_RNA_low),
    p_value = ifelse(
      n_pairs >= 3,
      t.test(
        x = GNG4_RNA_high,
        y = GNG4_RNA_low,
        paired = TRUE
      )$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_adj_BH)
scenic_tests_gng4

# Save results
write_xlsx(donor_scenic_df, path = '/filepath/fig5/fig5d/tonsil_gcliketfh_nfatc1_bach1_scenic_donor_scores.xlsx')
write_xlsx(scenic_tests_gng4, path = '/filepath/fig5/fig5d/tonsil_gcliketfh_nfatc1_bach1_scenic_donor_statistics.xlsx')

# Fig 5E - Experimental schematic for in vitro stimulation ----

# PDF exported from Biorender

# Fig 5F - Representative histograms of CD25, PD1, and Gγ4 protein expression in gated CD4 T cell subsets with versus without stimulation ----

# PDF exported from FlowJo, further annotated in Illustrator

# Fig 5G - Barplot of % Gγ4+ cells for each CD4 T cell subset in stimulated versus unstimulated conditions ----

# Import percentage values from FlowJo table export
tfh_gng4_flow_df <- read_excel("/filepath/fig5/fig5_flow_data.xlsx", sheet = "Gy4_Stim_Freq")

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
pdf('/filepath/fig5/fig5g/tfh_stim_gy4_freq_barplot.pdf', width = 5.75, height = 4)
tfh_stim_gy4_freq_barplot
dev.off()

# Statistics for % Gγ4+ differences between paired stimulated versus unstimulated conditions per donor were determined in Prism (Data file S11)
# Paired t-test with two-stage step-up Benjamini, Krieger, and Yekutieli FDR correction for multiple comparisons, *Q < 0.05, **Q < 0.01
# *Q and **Q symbols added in Illustrator

# Gγ4, PD1, and ICOS GMFI values for each subset were analyzed and visualized in Prism, shown in Fig. S13E-S13J

# Fig. 5H-5K - Flow cytometric analysis of Gγ4+ versus Gγ4- subsets of conventionally gated 'GC Tfh' ----

# Analyzed in FlowJo and Prism (Data file S11), visualization assembled in Illustrator