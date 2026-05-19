# Sam Barnett Dubensky et al.
# Derek A. Oldridge & Laura A. Vella Labs at the Children's Hospital of Philadelphia
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned Tfh in human lymphoid tissue
# Code and data visualization for Fig. 6
# Fig. 6 - GNG4 is induced in Tfh following vaccination and associated with pathogenic Tfh-like states in autoimmunity

# Set up R working environment ----

# Set working directory to Fig 6
setwd('/filepath/fig6/')

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
library(writexl) # 1.5.4
library(ggnewscale) # 0.5.0
library(readxl) # 1.4.3
library(UCell) # 2.8.0
library(scales) # 1.3.0
library(rlang) # 1.1.4
library(ggrepel) # 0.9.5
library(speckle) # 1.4.0

# Fig. 6A - Reanalysis of longitudinal SARS2 mRNA vaccination FNA scRNAseq study ----

# Reanalysis of Borcherding et al. CD4+ T cells exhibit distinct transcriptional phenotypes in the lymph nodes and blood following mRNA vaccination in humans. Nat Immunol. 2024. PMID 39164479. https://pubmed.ncbi.nlm.nih.gov/39164479/
# Import T cell subclustering Seurat object from SARS2 mRNA vaccination LN FNA scRNAseq study
# .rds file obtained from Figure2 tab of https://cellpilot.emed.wustl.edu/
covid_tcell_obj <- readRDS('/filepath/reanalysis/borcherding_mrna_vax/BorcherdingFig2.rds')

# Refer to Fig. S15 and related code in 'figS15_reanalysis_human_vaccination_scrnaseq_gng4_gctfh.R' for further annotation details of GC Tfh cluster

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
quantile(covid_gctfh_freq_df$percent_GCTfh, 0.25) - 1.5*IQR(covid_gctfh_freq_df$percent_GCTfh) # 4.274979 % GCTfh lower bound outlier threshold. 3.33% at d110 < 4.274979% lower bound threshold, suggesting that d110 may be an outlier

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
borcherding_dotplot_feats <- c('IL21','TOX2','BCL6','S1PR2','B3GAT1','NFATC1','BACH1','GNG4') # features of interest
pdf('/filepath/fig6/fig6a/borcherding_covid_fna_reanalysis_dotplot.pdf', width = 6.25, height = 6.5)
borcherding_covid_fna_reanalysis_dotplot <- DotPlot(covid_tcell_obj_filt, features = borcherding_dotplot_feats, group.by = 'cluster_timepoint', scale.by = 'size', dot.scale = 12, cols = 'RdBu') +
  theme(legend.justification.right = 'bottom',
        legend.title = element_text(size = 12, family = 'sans'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 14, face = 'italic', family = 'sans', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, family = 'sans')) +
  scale_y_discrete(labels = covid_tcell_labels) +
  scale_color_gradient2(
    low = '#4E79A7',           
    mid = "white",                
    high = 'coral3',        
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
borcherding_gctfh_time_bin_stats$data %>% filter(features.plot == 'GNG4') %>% select(features.plot,id,pct.exp,avg.exp) %>% View()

# d60 / d28 Pct Expr
36.3225806/27.2922636 # 1.330875 fold increase
# d201 / d60 Pct Expr
43.8006952/36.3225806 # 1.205881 fold increase

# d60 / d28 Avg Norm Expr
1.59654083/1.28730403 # 1.24022 fold increase
# d201 / d60 Avg Norm Expr
1.96361706/1.59654083 # 1.22992 fold increase

# Fig. 6B - Reanalysis of quadrivalent inactivated influenza longitudinal FNA scRNAseq data ----

# Reanalysis of Schattgen et al. Influenza vaccination stimulates maturation of the human T follicular helper cell response. Nat Immunol. 2024. PMID 39164477. https://pubmed.ncbi.nlm.nih.gov/39164477/
# Import Seurat objects for T cell clustering and Level 2 Tfh subclustering from IIV LN FNA scRNAseq study 
# Seurat object .rds files ('intergrated_Tcells_harmony_bydonor.rds' and 'intergrated_Tfh_harmony_bydonor.rds') obtained from Zenodo data repository (Version 3, Oct 20 2023, https://zenodo.org/records/12611325)
flu_tcell_obj <- readRDS('/filepath/reanalysis/intergrated_Tcells_harmony_bydonor.rds')
flu_tfh_obj <- readRDS('/filepath/reanalysis/intergrated_Tfh_harmony_bydonor.rds')

# Note - as provided, the Level 1 T cell object does not have a specifically annotated 'GC Tfh' cluster, whereas the Level 2 Tfh subclustering object does ('GC')
# Therefore, in steps below we examine expression of GNG4 and other GC-like Tfh features in both the L2 (Fig. S15D-F) and L1 objects (Fig. S15G-I), then map the 'GC' cluster from L2 to L1 (Fig. S15J-K)
# Refer to Figure S15 and related code in 'figS15_reanalysis_human_vaccination_scrnaseq_gng4_gctfh.R' for further annotation details of GC Tfh cluster

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
    low = '#4E79A7',
    mid = "white",
    high = 'coral3',
    midpoint = 0
  )

# Adjust legends for dotplot - further annotated in Illustrator
pdf('/filepath/fig6/fig6b/flu_fna_reanalysis_dotplot_smallerdots.pdf', width = 6.25, height = 6.5)
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
schattgen_gctfh_time_bin_stats$data %>% filter(features.plot == 'GNG4') %>% select(features.plot,id,pct.exp,avg.exp) %>% View()
# GNG4 Avg Expr d0-28 vs d60+ = 0.26100881 vs 1.17268623 (4.4928990328x)
# GNG4 % Expr d0-28 vs d60+ = 11.769416 vs 42.776524 (3.6345494118x)

# Fig. 6C - GNG4 ATAC coverage plot for tonsillar L4 GC vs nonGC-like Tfh groups, with rheumatoid arthritis fine-mapping GWAS variant positions annotated, GNG4 DAP regions highlighted ----

# Import L3 3WNN Tfh-like Seurat object from Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors)
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Annotate L3 Tfh-like cell clusters
tfh_subclust_names <- c(
  "0" = "Tfh-Circ",
  "1" = "Tcm", # Tcm will be excluded in Tfh analyses below, as well as any PBMC
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

# Filter Tcm and PBMC from L3 object to create tonsillar Tfh-only L4 object
l4_teaseq_tonsil_tfh_obj <- subset(l3_teaseq_tfh_obj, tfh_wnn_annot != 'Tcm' & hto.tissue == 'Tonsil')
ncol(l4_teaseq_tonsil_tfh_obj) # 2543 cells

# Finding all ATAC peaks within GNG4 region
DefaultAssay(l4_teaseq_tonsil_tfh_obj) <- 'ATAC'
tfh_annotation <- Annotation(l4_teaseq_tonsil_tfh_obj)
gng4_range <- tfh_annotation[tfh_annotation$gene_name == 'GNG4', ]
tfh_peak_ranges <- granges(l4_teaseq_tonsil_tfh_obj[["ATAC"]]) # get GRanges for all peaks within ATAC assay
tfh_peak_gng4 <- subsetByOverlaps(tfh_peak_ranges, gng4_range) # filter for peaks within GNG4
tfh_peak_gng4_df <- as.data.frame(tfh_peak_gng4)
tfh_peak_gng4_df # 14 total peak regions in GNG4 identified by cellranger-arc

# Find GNG4 DAP in tonsillar GC vs nonGC-like Tfh (no Tcm)
DefaultAssay(l4_teaseq_tonsil_tfh_obj) <- 'ATAC'
Idents(l4_teaseq_tonsil_tfh_obj) <- 'gc_vs_nongc_like'
gc_vs_nongc_tfh_dap <- FindMarkers(l4_teaseq_tonsil_tfh_obj, ident.1 = 'GC', ident.2 = 'nonGC', assay = 'ATAC', min.pct = 0, logfc.threshold = 0)
gc_vs_nongc_tfh_dap_gr <- StringToGRanges(rownames(gc_vs_nongc_tfh_dap), sep = c("-", "-"))
gc_vs_nongc_tfh_dap_closest_feat <- ClosestFeature(
  object = l4_teaseq_tonsil_tfh_obj,
  regions    = gc_vs_nongc_tfh_dap_gr,
  annotation = Annotation(l4_teaseq_tonsil_tfh_obj)
)
gc_vs_nongc_tfh_dap$closest_gene <- gc_vs_nongc_tfh_dap_closest_feat$gene_name
gc_vs_nongc_tfh_dap$transcript_id <- gc_vs_nongc_tfh_dap_closest_feat$tx_id
gc_vs_nongc_tfh_dap$atac_peak <- rownames(gc_vs_nongc_tfh_dap)
gng4_dap_gc_tfh <- gc_vs_nongc_tfh_dap %>% filter(closest_gene == 'GNG4') # 14 total peak regions in GNG4 identified by cellranger-arc
gng4_dap_gc_tfh %>% filter(p_val_adj < 0.05 & closest_gene == 'GNG4') # of these 14, only 3 were found to be significantly enriched in the GC-like L4 Tfh group (refer to Fig 5A analysis above)

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

# Set colors for cCREs of interest
# Red highlight of GNG4 dELS with increased accessibility in GC-like Tfh
# Darker grey for promoter-like sequence and Exons 1A + 1B region
# Lighter grey for other enhancer-like regions
rdbu_colors <- brewer.pal(11, "RdBu")
gng4_dap_gr$color <- c(rdbu_colors[3],'#4D4D4D','darkgrey','darkgrey','darkgrey')

# Having identified co-accessible and differentially accessible regions of GNG4 in L4 GC-like Tfh, next we considered whether any known variants from GWAS/eQTL studies are found near these regions of interest

# Ishigaki, Sakaue, Terao et al. Multi-ancestry genome-wide association analyses identify novel genetic mechanisms in rheumatoid arthritis. Nat. Genet. 2022. PMID 36333501. https://pubmed.ncbi.nlm.nih.gov/36333501/
# Fine-mapping GWAS results include set of variants mapped to GNG4 that are associated with RA diagnosis
# Lead variant annotated as rs1188620266 (chr1:235800357:CAA:C) in original study
# Due to changes in gnomAD versions since this study and our reanalysis, we have chosen to refer to this variant using the updated rs61512163 identifier (chr1:235637057:CAA:C), which we confirmed with the original study authors
# Refer to Methods text and Data File 10 for additional GWAS and variant details

# rsIDs for GNG4 variants identified in RA GWAS, where known
gng4_snps <- c(
  'rs61512163','rs12133886','rs35458456','rs34999365', # note rs61512163 rsID used here rather than rs1188620266
  'rs10926320','N/A','rs10802904','rs6429213',
  'rs12133526','N/A','N/A','rs7555242',
  'rs4391655'
)

# Positions of GNG4 variants in hg38, determined using Broad Institute Liftover Webtool (refer to Data File S10)
gng4_snp_pos <- c(
  235637057, 235638484, 235638486, 235638482,
  235638201, 235637653, 235643279, 235642088,
  235643355, 235635572, 235631575, 235635266,
  235637611
)

# rsID/positions were visualized in tonsillar L4 Tfh GNG4 locus coverage plot, then further annotated as below
# Within GNG4 dELS-containing DAP region - rs6429213 (235642088)
# Proximal to same GNG4 dELS-containing DAP region - rs10802904 (235643279) and rs12133526 (235643355)

# Assigning colors to each variant for exploration - each thin colored line was then used to guide placement of thicker dark red lines in Illustrator
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
  ranges   = IRanges(start = gng4_snp_pos - 10, end = gng4_snp_pos + 10), # widening variant positions to make colored lines visible to guide thicker line placement in Illustrator
  snp      = gng4_snps
)
gng4_snp_gr$color <- gng4_snp_gr_colors

# Combine GRanges objects for RA GWAS variants and GNG4 tonsillar L4 GC-like Tfh DAPs
combined_gr <- c(gng4_snp_gr,gng4_dap_gr)

# Determine visualization window for coverage plot
start_region <- min(start(combined_gr)) - 1000
end_region   <- max(start(combined_gr)) + 2000
plot_region <- paste0("chr1-", start_region, "-", end_region)

# Assemble coverage plot
l4_teaseq_tonsil_tfh_obj$gc_vs_nongc_like <- factor(l4_teaseq_tonsil_tfh_obj$gc_vs_nongc_like,  levels = c('GC','nonGC'))
gng4_snp_plot <- CoveragePlot(
  subset(l4_teaseq_tonsil_tfh_obj, subset = hto.tissue == 'Tonsil'),
  region = plot_region,
  links = FALSE,
  annotation = TRUE,
  peaks = FALSE,
  layer = "data",
  features = NULL,
  group.by = "gc_vs_nongc_like",
  region.highlight = combined_gr
) & 
  scale_fill_manual(values = c("GC" = "black", "nonGC" = "black"))

# Adjust axis titles and gene annotation labels
gng4_snp_plot[[1]][[1]] <- gng4_snp_plot[[1]][[1]] +
  theme(axis.ticks = element_blank()) +
  labs(y = NULL, x = NULL) # "Normalized Peak Signal (0-160)" y-axis label added in Illustrator

gng4_snp_plot[[1]][[2]] <- gng4_snp_plot[[1]][[2]] +
  scale_color_manual(values = c("black", "black")) +
  theme(axis.text.x = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Adjust track heights
gng4_snp_plot <- gng4_snp_plot +
  patchwork::plot_layout(heights = c(4, 0.6))

# Save plot
pdf('/filepath/fig6/fig6c/gng4_ra_gwas_coverageplot.pdf', width = 5, height = 3.75)
gng4_snp_plot
dev.off()

# Create genome annotation track for full GNG4 locus

# Get full GNG4 annotation
ann <- Annotation(l4_teaseq_tonsil_tfh_obj)
gng4_range <- ann[ann$gene_name == "GNG4"]

# Collapse ranges of transcript variants to select mimumum and maximum coordinates
gng4_full_gr <- GRanges(
  seqnames = as.character(unique(seqnames(gng4_range))[1]),
  ranges = IRanges(
    start = min(start(gng4_range)),
    end   = max(end(gng4_range))
  )
)

# Define coordinates of full GNG4 region with padding
gng4_full_start <- min(start(gng4_range)) - 10000
gng4_full_end   <- max(end(gng4_range)) + 10000

# Create character string of coordinates for complete GNG4 locus with padding
gng4_full_region_string <- paste0(
  as.character(seqnames(gng4_full_gr)),
  "-",
  gng4_full_start,
  "-",
  gng4_full_end
)

# Define coordinates of region spanning 13 GWAS variants with buffer to highlight in zoomed-out locus annotation
start_region <- min(start(combined_gr)) - 1000
end_region   <- max(start(combined_gr)) + 2000
zoom_highlight_gr <- GRanges(
  seqnames = seqnames(gng4_full_gr),
  ranges = IRanges(
    start = start_region,
    end   = end_region
  )
)
zoom_highlight_gr$color <- "#4D4D4D"

# Assemble genome annotation
DefaultAssay(l4_teaseq_tonsil_tfh_obj) <- "ATAC"
gng4_annotation_strip <- CoveragePlot(
  object = subset(l4_teaseq_tonsil_tfh_obj),
  group.by = 'gc_vs_nongc_like',
  region = gng4_full_region_string,
  links = FALSE,
  annotation = TRUE,
  peaks = FALSE,
  features = NULL,
  assay = "ATAC",
  region.highlight = zoom_highlight_gr # replace with 'gng4_snp_gr' created above to make zoomed-out GNG4 locus track with GWAS variants annotated
)
gng4_annotation_strip[[1]][[2]] <- gng4_annotation_strip[[1]][[2]] + scale_color_manual(values = c('black','black'))

gng4_annotation_strip
pdf('/filepath/fig6/fig6c/gng4_ra_gwas_annotation.pdf', width = 10, height = 6)
gng4_annotation_strip
dev.off()

# Data File S10 shows linkage disequilibrium analysis between these RA GWAS variants and the eQTL variants found in reanalysis of data from Schmiedel et al. Sci Immunol 2022.
# An additional variant of GNG4 (rs551907784) was reported in Verma et al. Science 2024 GWAS, but not found to be in LD with the eQTL variants, and was not visualized

# Fig. 6D - Recreate Al Souz et al. Seurat object ----

# Reanalysis of kidney tissue scRNAseq data from adult healthy donors versus systemic lupus erythematosus (SLE) patients 
# Originally published by Al Souz et al. (2026) Immunity
# https://pubmed.ncbi.nlm.nih.gov/42013863/
# Organ injury in systemic autoimmunity is mediated by stem-like CD8+ T cells arising from tissue-draining lymph nodes

# Retrieve scRNAseq files and metadata from Broad Single Cell Portal - https://singlecell.broadinstitute.org/single_cell/study/SCP3488/amp2-lupus-kidney-single-cell

# Define path to folders containing Al Souz et al. files
base_dir <- "/filepath/fig6/SCP3488/expression/69c1753e771a5b0296040b35/single_cell_portal/expr_raw"

# Assemble cell-by-gene count matrix
counts <- ReadMtx(
  mtx = file.path(base_dir, "matrix_raw.mtx.gz"),
  features = file.path(base_dir, "features_raw.tsv.gz"),
  cells = file.path(base_dir, "barcodes_raw.tsv.gz"),
  feature.column = 1
)

# Create Seurat object
kidney_sle_obj <- CreateSeuratObject(
  counts = counts,
  project = "SCP3488"
)

# Import metadata from Al Souz et al.
md <- fread("/filepath/fig6/SCP3488/metadata/single_cell_portal/metadata.txt.gz", data.table = FALSE)

# Match cell barcode rows between metadata and Seurat object
rownames(md) <- md$NAME
md <- md[colnames(kidney_sle_obj), , drop = FALSE]

# Add metadata to Seurat object
kidney_sle_obj <- AddMetaData(kidney_sle_obj, metadata = md)

# Import UMAP reduction coordinates from Al Souz et al.
umap_tnk <- fread("/filepath/fig6/SCP3488/cluster/single_cell_portal/umap_TNK.txt.gz", data.table = FALSE)
colnames(umap_tnk)
head(umap_tnk)

# Match UMAP dataframe cell barcodes to Seurat object
rownames(umap_tnk) <- umap_tnk$NAME
umap_df <- umap_tnk[colnames(kidney_sle_obj), c('X', 'Y'), drop = FALSE]
umap_mat <- as.matrix(sapply(umap_df, as.numeric))
rownames(umap_mat) <- colnames(kidney_sle_obj)
colnames(umap_mat) <- c("UMAPTNK_1", "UMAPTNK_2")

# Add UMAP reduction to Seurat object
kidney_sle_obj[["umap_tnk"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAPTNK_",
  assay = DefaultAssay(kidney_sle_obj)
)

# Inspect metadata
ncol(kidney_sle_obj) # 58555 cells
table(kidney_sle_obj$cellType) # includes myeloid celltypes - removed in steps below

# Exclude non-T/ILC/NK cell types from object
tnk_celltypes <- unique(as.character(kidney_sle_obj$cellType))
tnk_celltypes <- tnk_celltypes[grepl("^(T|NK)", tnk_celltypes)]
tnk_celltypes
kidney_sle_obj <- subset(kidney_sle_obj, subset = cellType %in% tnk_celltypes)

# Inspect metadata
ncol(kidney_sle_obj) # 34736 cells
table(kidney_sle_obj$cellType)

# Perform typical Seurat workflow of log-normalization, feature selection, and feature scaling
# Reflecting Al Souz et al. methods - "RNA data was normalized by LogNormalize (feature expression measurement divided by the total expression in each cell multiplied by a scale factor of 10,000 then natural log transformed) in Seurat v5. Top 2000 variable genes were determined using FindVariableFeatures with ‘vst’ selection method, and the normalized RNA data was scaled for variable genes using ScaleData."
kidney_sle_obj <- NormalizeData(kidney_sle_obj)
kidney_sle_obj <- FindVariableFeatures(kidney_sle_obj)
kidney_sle_obj <- ScaleData(kidney_sle_obj, features = rownames(kidney_sle_obj))

# Save recreated Seurat object
saveRDS(kidney_sle_obj, '/filepath/fig6/kidney_sle_seurat_obj_recreated.rds')

# Fig. 6D - UMAP of cell types in Al Souz et al. adult SLE kidney tissue scRNAseq study ----

# Recreate Seurat object using steps above or import if data already processed
kidney_sle_obj <- readRDS('/filepath/fig6/kidney_sle_seurat_obj_recreated.rds')

# Specify colors and labels for UMAP visualization
kidney_sle_obj$cellType <- gsub(
  "γδ",
  "gd",
  kidney_sle_obj$cellType
)
celltype_levels <- unique(as.character(kidney_sle_obj$cellType))
celltype_cols <- c(
  "#e3b420",
  "#657ab0",
  "#b6b7ba",
  "#dfcc78",
  "#728782",
  "#c9a3c1",
  "#bad4f4",
  "#479b73",
  "#F28E2B",
  "#e19dc6",
  "#c9744d",
  "#62adc3",
  "coral3", #ed7a77
  "#85a798",
  "#AEC7E8",
  "#c49980",
  "#DADAEB",
  "#8897d7",
  '#FFBB78',
  '#7F7F7F',
  '#e1d4c7',
  '#DBDB8D',
  '#B07AA1',
  '#59A14F',
  '#9E9AC8',
  '#C7E9C0'
)
celltype_cols <- setNames(celltype_cols[seq_along(celltype_levels)], celltype_levels)

# Assemble UMAP
pdf('/filepath/fig6/fig6d/kidney_sle_dimplot.pdf',
    height = 6,
    width = 9)
DimPlot(kidney_sle_obj, reduction = 'umap_tnk', group.by = 'cellType', label = TRUE, repel = TRUE, label.box = TRUE, label.size = 4, pt.size = 0.4, cols = celltype_cols) + 
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# Fig. 6E - UMAP with overlaid gene expression in kidney cell types ----

# Each FeaturePlot is made with and without legend, later assembled in Adobe Illustrator

# GNG4 without legend
pdf('/filepath/fig6/fig6e/kidney_sle_gng4_featureplot.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'GNG4',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','coral3')) + 
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# GNG4 with legend
pdf('/filepath/fig6/fig6e/kidney_sle_gng4_featureplot_withlegend.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'GNG4',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','coral3')) + 
  coord_fixed() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# CXCL13 without legend
pdf('/filepath/fig6/fig6e/kidney_sle_cxcl13_featureplot.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'CXCL13',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','#4E79A7')) +
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# CXCL13 with legend
pdf('/filepath/fig6/fig6e/kidney_sle_cxcl13_featureplot_withlegend.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'CXCL13',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','#4E79A7')) + 
  coord_fixed() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# IL21 without legend
pdf('/filepath/fig6/fig6e/kidney_sle_il21_featureplot.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'IL21',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','#4E79A7')) + 
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# IL21 with legend
pdf('/filepath/fig6/fig6e/kidney_sle_il21_featureplot_withlegend.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'IL21',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','#4E79A7')) + 
  coord_fixed() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# CXCR5 without legend
pdf('/filepath/fig6/fig6e/kidney_sle_cxcr5_featureplot.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'CXCR5',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','#4E79A7')) + 
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# CXCR5 with legend
pdf('/filepath/fig6/fig6e/kidney_sle_cxcr5_featureplot_withlegend.pdf',
    height = 4,
    width = 6)
FeaturePlot(kidney_sle_obj, features = 'CXCR5',  reduction = 'umap_tnk', label = FALSE, pt.size = 1, order = TRUE, cols = c('#e2e2e2','#4E79A7')) + 
  coord_fixed() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# Fig. 6F - Comparing % GNG4+ cells of all cells per donor in healthy versus SLE kidney ----

# Define GNG4 RNA+ cells
kidney_sle_obj$GNG4_pos <- FetchData(kidney_sle_obj, vars = "GNG4")[, 1] > 0

# Get Seurat object metadata
kidney_sle_md <- kidney_sle_obj@meta.data

# Compute percent GNG4 RNA+ cells per donor out of all cells in kidney_sle_obj
donor_gng4_df <- kidney_sle_md %>%
  dplyr::group_by(
    donor_id,
    Type
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_gng4_pos = sum(GNG4_pos == "TRUE", na.rm = TRUE),
    pct_gng4_pos = 100 * n_gng4_pos / n_total,
    .groups = "drop"
  )

# Assemble plot
pdf(
  '/filepath/fig6/fig6f/kidney_sle_pct_gng4_pos_anycell_of_tnk_vlnplot.pdf',
  height = 3,
  width = 5
)
ggplot(
  donor_gng4_df,
  aes(x = Type, y = pct_gng4_pos, fill = Type)) +
  geom_violin() +
  geom_jitter(
    width = 0.25,
    height = 0,
    size = 5,
    shape = 21,
    color = "black",
    stroke = 0.5) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "black"
  ) +
  theme_classic() +
  labs(
    y = expression("% " * italic("GNG4") * "+ cells of all T/NK/ILC")
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "grey",
      "LN" = "coral3"
    )
  )
dev.off()

# Fig. 6F - Comparing % GNG4+ Tfh/Tph-like cells of all cells per donor in healthy versus SLE kidney ----

# Compute percent GNG4 RNA+ T14. Tfh/Tph cells per donor out of all cells in kidney_sle_obj
donor_gng4_df <- kidney_sle_md %>%
  dplyr::group_by(
    donor_id,
    Type
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_gng4_pos = sum(cellType == "T14. Tfh/Tph" & GNG4_pos == "TRUE", na.rm = TRUE),
    pct_gng4_pos = 100 * n_gng4_pos / n_total,
    .groups = "drop"
  )

# Plot
pdf('/filepath/fig6/fig6f/kidney_sle_pct_gng4_pos_tfh_of_tnk_vlnplot.pdf',
    height = 3,
    width = 5)
ggplot(
  donor_gng4_df,
  aes(x = Type, y = pct_gng4_pos, fill = Type))  +
  geom_violin() +
  geom_jitter(
    width = 0.25,
    size = 5,
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "black"
  ) +
  theme_classic() +
  labs(
    y = expression("% " * italic("GNG4") * "+ Tfh/Tph of all T/NK/ILC")
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "grey",
      "LN" = "coral3"
    )
  )
dev.off()

# Fig. 6F - Propeller tests for % GNG4+ cells and % GNG4+ Tfh/Tph-like cells of all cells in healthy versus SLE kidney tissue groups ----

# Prepare metadata for Propeller test
propeller_md <- kidney_sle_md %>%
  dplyr::filter(Type %in% c("Control", "LN")) %>% # Comparing healthy donors versus lupus nephritis (LN)
  dplyr::mutate(
    Type = factor(Type, levels = c("Control", "LN")), 
    
    # Metadata for Propeller test 1 (% GNG4+ any cell type of all cells per donor in healthy donors versus LN)
    # Classify each cell as GNG4+ or other
    gng4_anycell_class = dplyr::case_when(
      GNG4_pos %in% c(TRUE, "TRUE") ~ "GNG4_pos",
      TRUE ~ "Other_TNK"
    ),
    
    # Metadata for Propeller test 2 (% GNG4+ Tfh/Tph-like cell of all cells per donor in healthy donors versus LN)
    # Classify each cell as GNG4+ Tfh/Tph or other
    gng4_tfh_class = dplyr::case_when(
      cellType == "T14. Tfh/Tph" &
        GNG4_pos %in% c(TRUE, "TRUE") ~ "GNG4_pos_Tfh",
      TRUE ~ "Other_TNK"
    )
  )

# Check donor totals per group
propeller_md %>%
  dplyr::distinct(donor_id, Type) %>%
  dplyr::count(Type)

# Perform Propeller test 1 (% GNG4+ any cell type of all cells per donor in Control versus LN groups)
gng4_anycell_propeller <- propeller(
  clusters = propeller_md$gng4_anycell_class,
  sample = propeller_md$donor_id,
  group = propeller_md$Type,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)
gng4_anycell_propeller
gng4_anycell_propeller["GNG4_pos", , drop = FALSE]

# Perform Propeller test 2 (% GNG4+ Tfh/Tph of all cells per donor in Control versus LN groups)
gng4_tfh_propeller <- propeller(
  clusters = propeller_md$gng4_tfh_class,
  sample = propeller_md$donor_id,
  group = propeller_md$Type,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)
gng4_tfh_propeller
gng4_tfh_propeller["GNG4_pos_Tfh", , drop = FALSE]

# Fig. 6G - UMAP of cell types in Balasubramanian et al. pediatric SLE PBMC scRNAseq study ----

# Reanalysis of Balasubramanian et al. (2025) Nat Immunol
# Single-cell RNA profiling of blood CD4+ T cells identifies distinct helper and dysfunctional regulatory clusters in children with SLE
# https://pubmed.ncbi.nlm.nih.gov/41120754/
# Retrieved 'sorted.cd4.total.rds' file from Zenodo repository - https://zenodo.org/records/16385794

# Import Seurat object
pbmc_sle_obj <- readRDS('/filepath/fig6/sorted.cd4.total.rds')

# Annotation of activated cTfh-like cluster (based on reanalysis including features in Fig. 6H)
pbmc_sle_obj$clusters_grouped <- gsub(
  "Proliferating",
  "Act. cTfh-like",
  pbmc_sle_obj$clusters_grouped
)

# Define GNG4+ cells
pbmc_sle_obj$GNG4_pos <- FetchData(pbmc_sle_obj, vars = "GNG4")[,1] > 0

# Set colors for clusters
celltype_levels <- unique(as.character(pbmc_sle_obj$clusters_grouped))
celltype_cols <- c(
  "#657ab0",
  "#85a798",
  "#c9a3c1",
  "coral3",
  "#dfcc78"
)
celltype_cols <- setNames(celltype_cols[seq_along(celltype_levels)], celltype_levels)

# UMAP with cluster labels
pdf('/filepath/fig6/fig6g/pediatric_pbmc_sle_dimplot.pdf',
    height = 5,
    width = 7)
DimPlot(pbmc_sle_obj, group.by = 'clusters_grouped', raster = FALSE, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 5, pt.size = 0.2, cols = celltype_cols) + 
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# UMAP without labels, applied in Illustrator
pdf('/filepath/fig6/fig6g/pediatric_pbmc_sle_dimplot_nolabs.pdf',
    height = 5,
    width = 7)
DimPlot(pbmc_sle_obj, group.by = 'clusters_grouped', raster = FALSE, pt.size = 0.2, cols = celltype_cols) + 
  coord_fixed() + NoLegend() + theme(plot.title = element_blank()) + NoAxes()
dev.off()

# Fig. 6H - Dotplot of Tfh-associated genes across clusters in SLE PBMC study ----

# Specify order of clusters
pbmc_sle_obj$clusters_grouped <- factor(
  pbmc_sle_obj$clusters_grouped,
  levels = c(
    "Naive",
    "ISG-high",
    "Treg",
    "Memory",
    "Act. cTfh-like"
  )
)
Idents(pbmc_sle_obj) <- 'clusters_grouped'

# Assemble dotplot
pdf('/filepath/fig6/fig6h/pbmc_sle_gng4_annotation_dotplot.pdf',
    height = 4,
    width = 6.5)
DotPlot(pbmc_sle_obj, features = c('CXCR5','B3GAT1','BTLA','CXCL13','IL21','IL4','BCL6','TOX2','POU2AF1','ASCL2','NFATC1','GNG4'), 
        group.by = 'clusters_grouped', dot.scale = 16, cluster.idents = FALSE) +
  scale_color_gradient2(low = "#4E79A7",mid = "white",high = "coral3",midpoint = 0) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
dev.off()

# Fig. 6I - Comparing % GNG4+ cells of all cells per donor in PBMC from healthy vs SLE w/o LN vs SLE w/ LN donors ----

# Define GNG4 RNA+ cells
pbmc_sle_obj$GNG4_pos <- FetchData(pbmc_sle_obj, vars = "GNG4")[,1] > 0

# Extract metadata
pbmc_sle_md <- pbmc_sle_obj@meta.data

# Order groups for visualization
pbmc_sle_md$Group_LN <- factor(
  pbmc_sle_md$Group_LN,
  levels = c("HD", "NoLN", "LN")
)

# Determine % GNG4 RNA+ cells of all CD4 T cells per donor by group
donor_gng4_anycell_df <- pbmc_sle_md %>%
  dplyr::group_by(
    Subject_ID,
    Group_LN
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_gng4_pos = sum(
      GNG4_pos %in% c(TRUE, "TRUE"),
      na.rm = TRUE
    ),
    pct_gng4_pos = 100 * n_gng4_pos / n_total,
    .groups = "drop"
  )

# Assemble violin plot
pdf(
  "/filepath/fig6/fig6i/sle_cd4_percent_gng4_pos_Group_LN_vlnplot.pdf",
  height = 2.5,
  width = 5
)

ggplot(
  donor_gng4_anycell_df,
  aes(x = Group_LN, y = pct_gng4_pos, fill = Group_LN)
) +
  geom_violin() +
  geom_jitter(
    width = 0.25,
    size = 5,
    height = 0,
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "black"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "HD" = "grey",
      "NoLN" = "coral",
      "LN" = "coral3"
    )
  )

dev.off()

# Fig. 6I - Comparing % GNG4+ activated Tfh-like cells of all cells per donor in PBMC from healthy vs SLE w/o LN vs SLE w/ LN donors ----

# Determine % GNG4 RNA+ activated cTfh-like cells of all CD4 T cells per donor by group
donor_gng4_tfh_df <- pbmc_sle_md %>%
  dplyr::group_by(
    Subject_ID,
    Group_LN
  ) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_gng4_pos = sum(
      clusters_grouped == "Act. cTfh-like" & GNG4_pos %in% c(TRUE, "TRUE"),
      na.rm = TRUE
    ),
    pct_gng4_pos = 100 * n_gng4_pos / n_total,
    .groups = "drop"
  )

# Assemble violin plot
pdf(
  "/filepath/fig6/fig6i/pbmc_sle_percent_gng4_pos_act_ctfh_Group_LN_vlnplot.pdf",
  height = 2.5,
  width = 5
)

ggplot(
  donor_gng4_tfh_df,
  aes(x = Group_LN, y = pct_gng4_pos, fill = Group_LN)
) +
  geom_violin() +
  geom_jitter(
    width = 0.25,
    height = 0,
    size = 5,
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "black"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "HD" = "grey",
      "NoLN" = "coral1",
      "LN" = "coral3"
    )
  )

dev.off()

# Fig. 6I - Pairwise Propeller tests for % GNG4+ CD4 T cells and % GNG4+ activated cTfh-like cells of all CD4 T cells in PBMC from healthy vs SLE w/o LN vs SLE w/ LN donors ----

# Prepare metadata for pairwise Propeller tests
propeller_md <- pbmc_sle_md %>%
  dplyr::filter(Group_LN %in% c("HD", "NoLN", "LN")) %>% # Comparing healthy donors versus SLE w/o LN versus SLE w/ LN
  dplyr::mutate(
    Group_LN = factor(Group_LN, levels = c("HD", "NoLN", "LN")),
    
    # Metadata for pairwise Propeller test 1 (% GNG4+ cells of all CD4 T cells per donor in HD versus SLE w/o LN versus SLE w/ LN)
    # Classify each cell as GNG4+ CD4 T cell or other CD4 T cell
    gng4_anycell_class = dplyr::case_when(
      GNG4_pos %in% c(TRUE, "TRUE") ~ "GNG4_pos",
      TRUE ~ "Other_CD4"
    ),
    
    # Metadata for pairwise Propeller test 2 (% GNG4+ activated cTfh-like cells of all CD4 T cells per donor in HD versus SLE w/o LN versus SLE w/ LN)
    # Classify each cell as GNG4+ activated cTfh-like or other CD4 T cell
    gng4_act_ctfh_class = dplyr::case_when(
      clusters_grouped == "Act. cTfh-like" &
        GNG4_pos %in% c(TRUE, "TRUE") ~ "GNG4_pos_Act_cTfh",
      TRUE ~ "Other_CD4"
    )
  )

# Inspect donor totals per group
propeller_md %>%
  dplyr::distinct(Subject_ID, Group_LN) %>%
  dplyr::count(Group_LN)

# Perform pairwise Propeller test 1 (% GNG4+ any CD4 T cell type of all CD4 T cells per donor in HD vs SLE w/o LN vs SLE w/ LN)

# HD vs SLE w/o LN
propeller_md_hd_noln <- propeller_md %>%
  dplyr::filter(Group_LN %in% c("HD", "NoLN")) %>%
  dplyr::mutate(
    Group_LN = factor(Group_LN, levels = c("HD", "NoLN"))
  )

gng4_anycell_propeller_hd_noln <- propeller(
  clusters = propeller_md_hd_noln$gng4_anycell_class,
  sample = propeller_md_hd_noln$Subject_ID,
  group = propeller_md_hd_noln$Group_LN,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)

# HD vs SLE w/ LN
propeller_md_hd_ln <- propeller_md %>%
  dplyr::filter(Group_LN %in% c("HD", "LN")) %>%
  dplyr::mutate(
    Group_LN = factor(Group_LN, levels = c("HD", "LN"))
  )

gng4_anycell_propeller_hd_ln <- propeller(
  clusters = propeller_md_hd_ln$gng4_anycell_class,
  sample = propeller_md_hd_ln$Subject_ID,
  group = propeller_md_hd_ln$Group_LN,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)

# SLE w/o LN vs SLE w/ LN
propeller_md_noln_ln <- propeller_md %>%
  dplyr::filter(Group_LN %in% c("NoLN", "LN")) %>%
  dplyr::mutate(
    Group_LN = factor(Group_LN, levels = c("NoLN", "LN"))
  )

gng4_anycell_propeller_noln_ln <- propeller(
  clusters = propeller_md_noln_ln$gng4_anycell_class,
  sample = propeller_md_noln_ln$Subject_ID,
  group = propeller_md_noln_ln$Group_LN,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)

# Extract GNG4+ CD4 T cell row from each pairwise Propeller test and perform BH-adjustment across pairwise tests
gng4_anycell_pairwise_propeller <- dplyr::bind_rows(
  data.frame(
    comparison = "HD vs NoLN",
    gng4_anycell_propeller_hd_noln["GNG4_pos", , drop = FALSE]
  ),
  data.frame(
    comparison = "HD vs LN",
    gng4_anycell_propeller_hd_ln["GNG4_pos", , drop = FALSE]
  ),
  data.frame(
    comparison = "NoLN vs LN",
    gng4_anycell_propeller_noln_ln["GNG4_pos", , drop = FALSE]
  )
) %>%
  dplyr::mutate(
    p_adj_BH_across_pairwise = p.adjust(P.Value, method = "BH")
  )

gng4_anycell_pairwise_propeller

# Perform pairwise Propeller test 2 (% GNG4+ activated cTfh-like cells of all CD4 T cells per donor in HD vs SLE-LN vs SLE+LN)

# HD vs SLE w/o LN
gng4_act_ctfh_propeller_hd_noln <- propeller(
  clusters = propeller_md_hd_noln$gng4_act_ctfh_class,
  sample = propeller_md_hd_noln$Subject_ID,
  group = propeller_md_hd_noln$Group_LN,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)

# HD vs SLE w/ LN
gng4_act_ctfh_propeller_hd_ln <- propeller(
  clusters = propeller_md_hd_ln$gng4_act_ctfh_class,
  sample = propeller_md_hd_ln$Subject_ID,
  group = propeller_md_hd_ln$Group_LN,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)

# SLE w/o LN vs SLE w/ LN
gng4_act_ctfh_propeller_noln_ln <- propeller(
  clusters = propeller_md_noln_ln$gng4_act_ctfh_class,
  sample = propeller_md_noln_ln$Subject_ID,
  group = propeller_md_noln_ln$Group_LN,
  transform = "logit",
  robust = TRUE,
  trend = FALSE
)

# Extract GNG4+ activated cTfh-like row from each pairwise Propeller test and perform BH-adjustment across tests
gng4_act_ctfh_pairwise_propeller <- dplyr::bind_rows(
  data.frame(
    comparison = "HD vs NoLN",
    gng4_act_ctfh_propeller_hd_noln["GNG4_pos_Act_cTfh", , drop = FALSE]
  ),
  data.frame(
    comparison = "HD vs LN",
    gng4_act_ctfh_propeller_hd_ln["GNG4_pos_Act_cTfh", , drop = FALSE]
  ),
  data.frame(
    comparison = "NoLN vs LN",
    gng4_act_ctfh_propeller_noln_ln["GNG4_pos_Act_cTfh", , drop = FALSE]
  )
) %>%
  dplyr::mutate(
    p_adj_BH_across_pairwise = p.adjust(P.Value, method = "BH")
  )
gng4_act_ctfh_pairwise_propeller