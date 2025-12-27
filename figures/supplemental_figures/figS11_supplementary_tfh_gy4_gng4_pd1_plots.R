# Barnett Dubensky et al. 2025 bioRxiv
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned CD4 T follicular helper cells in humans
# Code and data visualization for relevant panels of Fig. S11
# Fig. S11 - Gating strategy and quantification of Gγ4 protein and GNG4 RNA expression in human CD4 T cells. 

# Set up working environment ----

# Set working directory
setwd('/filepath/fig3/fig3_supp/')

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
library(readxl) # 1.4.3

# Specify colors for each tissue
tissue_cols <- c("PBMC" = "#4b94d6", "Tonsil" = "#e49a78")

# Import Seurat object from TEAseq Data Preprocessing Step 15 (trimodal dimensionality reduction and L3 3WNN subclustering of Tfh-like cells from L2 T cell object, including Harmony integration across donors) ----
l3_teaseq_tfh_obj <- readRDS('/filepath/step15_tfh_subcluster/l3_teaseq_tfh_obj.rds')

# Filter Tcm from L3 object to create Tfh-only L4 object ----
l4_teaseq_tfh_no_tcm_obj <- subset(l3_teaseq_tfh_obj, tfh_wnn_annot != 'Tcm')

# Fig S11C-S11D, supplement to Fig. 3J -  Gy4 GMFI in Tfh subsets vs nonTfh and Naive cells ----

# Tfh Subsets Gy4 GMFI Barplot
tfh_gng4_flow_df <- read_excel("/filepath/fig3/fig3_flow_data.xlsx", sheet = "Gy4_GMFI")

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
    values_to = "GMFI"                     
  )

tfh_gng4_flow_df_long <- tfh_gng4_flow_df_long %>%
  mutate(
    Subset = factor(
      Subset,
      levels = c("Naïve", "nonTfh", "PD1-", "PD1+", "PD1br")
    )
  )

# Tonsil sample Gy4 GMFI barplot with 95% confidence intervals
tfh_gng4_flow_df_ton <- tfh_gng4_flow_df_long %>% subset(Tissue == 'Tonsil')
ton_tfh_gy4_gmfi_barplot <- ggplot(tfh_gng4_flow_df_ton, aes(x = Subset, y = GMFI, fill = Subset)) +
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
  scale_y_continuous(breaks = seq(0, 5500, 500)) +
  coord_cartesian(ylim = c(0, 5500))
ton_tfh_gy4_gmfi_barplot
pdf('/filepath/fig3/fig3_supp/ton_tfh_gy4_gmfi_barplot.pdf', width = 5, height = 5)
ton_tfh_gy4_gmfi_barplot
dev.off()

# PBMC sample Gy4 GMFI barplot with 95% confidence intervals
tfh_gng4_flow_df_pbmc <- tfh_gng4_flow_df_long %>% subset(Tissue == 'PBMC')
pbmc_tfh_gy4_gmfi_barplot <- ggplot(tfh_gng4_flow_df_pbmc, aes(x = Subset, y = GMFI, fill = Subset)) +
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
  scale_y_continuous(breaks = seq(0, 5500, 500)) +
  coord_cartesian(ylim = c(0, 5500))
pbmc_tfh_gy4_gmfi_barplot
pdf('/filepath/fig3/fig3_supp/pbmc_tfh_gy4_gmfi_barplot.pdf', width = 5, height = 5)
pbmc_tfh_gy4_gmfi_barplot
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
  path = "/filepath/fig3/fig3_supp/tfh_gy4_freq_gmfi_95ci_stats.xlsx"
)

all_stats_t <- as.data.frame(t(all_stats))
all_stats_t$variables <- rownames(all_stats_t)
all_stats_t <- all_stats_t %>% select(variables, everything())
writexl::write_xlsx(
  all_stats_t,
  path = "/filepath/fig3/fig3_supp/tfh_gy4_freq_gmfi_95ci_stats_transpose.xlsx"
)

# Fig S11G - Finding % PD1br Tfh of all Tfh in Tonsil or PBMC ----

# RidgePlot visualization of PD1 expression in tonsil versus PBMC Tfh
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'hto.tissue'
l4_teaseq_tfh_no_tcm_obj$hto.tissue <- fct_relevel(l4_teaseq_tfh_no_tcm_obj$hto.tissue, "Tonsil", "PBMC")
pdf('/filepath/fig3/fig3_supp/pd1_ridgeplot_tfh.pdf', width = 3.5, height = 3.5)
RidgePlot(l4_teaseq_tfh_no_tcm_obj, features = 'adt-CD279', group.by = 'hto.tissue') + 
  scale_fill_manual(values = tissue_cols) + 
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank()) +
  NoLegend() +
  geom_vline(xintercept = 1.7) +
  geom_vline(xintercept = 2.35)
dev.off()

# Finding number of Tfh in each PD1 bin for PBMC versus tonsil
l4_tfh_md <- l4_teaseq_tfh_no_tcm_obj@meta.data
l4_teaseq_tfh_no_tcm_obj <- JoinLayers(l4_teaseq_tfh_no_tcm_obj, assay = 'ADT')
tfh_pd1_vals <- GetAssayData(l4_teaseq_tfh_no_tcm_obj, assay = "ADT", slot = "data")["adt-CD279", colnames(l4_teaseq_tfh_no_tcm_obj)]
tfh_pd1_df <- tibble(
  Cell        = colnames(l4_teaseq_tfh_no_tcm_obj),
  Tissue      = l4_teaseq_tfh_no_tcm_obj$hto.tissue,
  PD1   = as.numeric(tfh_pd1_vals)
) %>%
  mutate(Bin = dplyr::case_when(
    PD1 < 1.7            ~ "PD1-",
    PD1 >= 1.7 & PD1 < 2.35 ~ "PD1+",
    PD1 >= 2.35          ~ "PD1br",
    TRUE ~ NA_character_
  ))
tfh_pd1_df
table(tfh_pd1_df$Bin, tfh_pd1_df$Tissue) 
# PBMC - 3 PD1br / 3+81+369 total = 0.66% PD1br of all Tfh across PBMC samples
# Tonsil - 1707 PD1br / 1707+421+415 total = 67.1% PD1br of all Tfh across tonsil samples

# Finding percentage of PD1br Tfh per donor in each tissue - visualized by scatter plot in Prism
l4_teaseq_tfh_no_tcm_obj$pd1_adt_class <- ifelse(FetchData(l4_teaseq_tfh_no_tcm_obj, vars = "adt-CD279") >= 2.35, "PD1_ADT_br", "PD1_ADT_other")
pct_high_by_donor <- l4_teaseq_tfh_no_tcm_obj@meta.data %>%
  transmute(
    donor = hto.donor,
    class = trimws(as.character(pd1_adt_class))
  ) %>%
  filter(!is.na(donor), !is.na(class)) %>%
  group_by(donor) %>%
  summarize(
    n_total = n(),
    n_high  = sum(class == "PD1_ADT_br"),
    pct_high = 100 * n_high / n_total,
    .groups = "drop"
  ) %>%
  arrange(donor)
pct_high_by_donor

# Fig S11H - Finding % GNG4 RNA+ Tfh of all Tfh in Tonsil or PBMC ----

# RidgePlot visualization of GNG4 RNA expression in tonsil versus PBMC Tfh
Idents(l4_teaseq_tfh_no_tcm_obj) <- 'hto.tissue'
DefaultAssay(l4_teaseq_tfh_no_tcm_obj) <- 'RNA'
l4_teaseq_tfh_no_tcm_obj$hto.tissue <- fct_relevel(l4_teaseq_tfh_no_tcm_obj$hto.tissue, "PBMC", "Tonsil")
pdf('/filepath/fig3/fig3_supp/gng4_vlnplot_tfh.pdf', width = 3.5, height = 3.5)
VlnPlot(l4_teaseq_tfh_no_tcm_obj, features = 'rna_GNG4', group.by = 'hto.tissue', layer = 'data') + 
  scale_fill_manual(values = tissue_cols) + 
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank()) +
  NoLegend() +
  geom_hline(yintercept = 1) # specified threshold for RNA+ versus RNA-
dev.off()

# Finding number of Tfh in each GNG4 RNA bin for PBMC vs Tonsil
l4_teaseq_tfh_no_tcm_obj$gng4_rna_class <- ifelse(FetchData(l4_teaseq_tfh_no_tcm_obj, vars = "rna_GNG4") >= 1, "GNG4_RNA_high", "GNG4_RNA_low")
table(l4_teaseq_tfh_no_tcm_obj$gng4_rna_class, l4_teaseq_tfh_no_tcm_obj$hto.tissue) 
# 9/(9+444)*100 PBMC = 1.9867549669% of PBMC cTfh GNG4 RNA+
# 588/(588+1955)*100 Tonsil = 23.1222965002% of Tonsil Tfh are GNG4 RNA+

# Finding percentage of GNG4 RNA+ Tfh per donor in each tissue - visualized by scatter plot in Prism
pct_high_by_donor <- l4_teaseq_tfh_no_tcm_obj@meta.data %>%
  transmute(
    donor = hto.donor,
    class = trimws(as.character(gng4_rna_class))
  ) %>%
  filter(!is.na(donor), !is.na(class)) %>%
  group_by(donor) %>%
  summarize(
    n_total = n(),
    n_high  = sum(class == "GNG4_RNA_high"),
    pct_high = 100 * n_high / n_total,
    .groups = "drop"
  ) %>%
  arrange(donor)
pct_high_by_donor

# Remainder of plots in Fig S11 derive from FlowJo and Prism visualization