# Step 10 - Quality Control & Filtering (GEM5) ----

# Substeps
# 1) Set up R working environment
# 2) Import Seurat object for GEM well (empty-filtered, HTO-demultiplexed, SNP-demultiplexed, and multiplet-filtered)
# 3) Calculate RNA QC metrics
# 4) Set RNA QC thresholds
# 5) Calculate ATAC QC metrics
# 6) Set ATAC QC thresholds
# 7) Calculate ADT QC metrics
# 8) Select ADT QC thresholds
# 9) Calculate HTO QC metrics
# 10) Select HTO QC thresholds
# 11) Save Seurat object with all QC metrics appended before filtering
# 12) Filter Seurat object based on QC thresholds selected above for each modality
# 13) Assess number of cells filtered per sample
# 14) Save new object after QC filtering
# 15) Repeat for all GEM wells and proceed to next step

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
library(clustree) # 0.5.1
library(SoupX) # 1.6.2
library(scDblFinder) # 1.19.1

# Set working directory for this step
setwd("/filepath/step10_quality_control")

# 2) Import Seurat object for GEM well (empty-filtered, HTO-demultiplexed, SNP-demultiplexed, and multiplet-filtered) ----
gem5.raw <- readRDS("/filepath/step9_annotate_snp/gem5.annot.snp.rds") 

# Import quality control metric csv from cellranger-arc output for each GEM well
gem5.arc.metrics <- read.csv('/filepath/cellranger-arc_output/GEM5_ARC/outs/per_barcode_metrics.csv')

# Add GEM5 metadata
rownames(gem5.arc.metrics) <- gem5.arc.metrics$barcode
colnames(gem5.arc.metrics)[colnames(gem5.arc.metrics) == "barcode"] <- "arc_barcode"
gem5.rownames <- intersect(rownames(gem5.arc.metrics), colnames(gem5.raw))
gem5.arc.md <- gem5.arc.metrics[gem5.rownames,]
gem5.raw <- AddMetaData(gem5.raw, gem5.arc.md)

# 3) Calculate RNA QC Metrics ----

DefaultAssay(gem5.raw) <- 'RNA'
summary(gem5.raw$nCount_RNA) # RNA Count median -
summary(gem5.raw$nFeature_RNA) # RNA Feature median - 
summary(gem5.raw$percent.mt) # % Mitochondrial gene UMI median - 

# Calculate % Ribosomal gene UMI
str_view(rownames(gem5.raw), pattern = "^RP[LS]") # check gene names with prefix RP[LS]
gem5.raw[["percent.ribo"]] <- PercentageFeatureSet(gem5.raw, pattern = "^RP[LS]") # calculate % of UMI from ribosomal genes
summary(gem5.raw$percent.ribo) # % Ribosomal gene UMI median - 

# Calculate % Hemoglobin gene UMI
hbgene.list <- c('HBB','HBD','HBG1','HBG2','HBE1','HBZ','HBM','HBA2','HBA1','HBQ1') # curated list of hemoglobin-related genes
gem5.raw[["percent.hb"]] <- PercentageFeatureSet(gem5.raw, features = hbgene.list) # calculate % of UMI from hemoglobin-related genes
summary(gem5.raw$percent.hb) # % Hemoglobin gene UMI median - 

# Calculate % Platelet gene UMI
gem5.raw[["percent.plt"]] <- PercentageFeatureSet(gem5.raw, pattern = "PF4|PPBP") # calculate % of UMI from platelet-related genes (select few, relatively specific)
summary(gem5.raw$percent.plt) # % Platelet gene UMI median -

# 4) Select RNA QC thresholds ----

# Compare QC threshold across sorting groups
Idents(gem5.raw) <- "hto.sort"

# RNA Count
VlnPlot(gem5.raw, features = "nCount_RNA", pt.size = 0.1, log = TRUE) + NoLegend() & geom_hline(yintercept = c(600,12500))
VlnPlot(gem5.raw, features = "nCount_RNA", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(600,12500))
RidgePlot(gem5.raw, features = "nCount_RNA", log = TRUE) + NoLegend() & geom_vline(xintercept = c(600,12500))
RidgePlot(gem5.raw, features = "nCount_RNA", log = FALSE) + NoLegend() & geom_vline(xintercept = c(600,12500))

# RNA Feature
VlnPlot(gem5.raw, features = "nFeature_RNA", pt.size = 0.1, log = TRUE) + NoLegend() & geom_hline(yintercept = c(450,4000))
VlnPlot(gem5.raw, features = "nFeature_RNA", pt.size = 0.1, log = FALSE, y.max = 5100) + NoLegend() & geom_hline(yintercept = c(450,4000))
RidgePlot(gem5.raw, features = "nFeature_RNA", log = TRUE) + NoLegend() & geom_vline(xintercept = c(450,4000))
RidgePlot(gem5.raw, features = "nFeature_RNA", log = FALSE) + NoLegend() & geom_vline(xintercept = c(450,4000))

# RNA Count vs Feature
FeatureScatter(gem5.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", log = TRUE) & 
  geom_hline(yintercept = c(450,4000)) & 
  geom_vline(xintercept = c(600,12500))
DensityScatter(gem5.raw, x = "nCount_RNA", y = "nFeature_RNA", log_x = TRUE, log_y = TRUE) & 
  geom_hline(yintercept = c(450,4000)) & 
  geom_vline(xintercept = c(600,12500))

# % Mitochondrial UMI
VlnPlot(gem5.raw, features = "percent.mt", pt.size = 0.1, log = TRUE) + NoLegend() & geom_hline(yintercept = c(25))
VlnPlot(gem5.raw, features = "percent.mt", pt.size = 0.1, log = FALSE, y.max = 31) + NoLegend() & geom_hline(yintercept = c(25))
RidgePlot(gem5.raw, features = "percent.mt", log = TRUE) + NoLegend() & geom_vline(xintercept = c(25))
RidgePlot(gem5.raw, features = "percent.mt", log = FALSE) + NoLegend() & geom_vline(xintercept = c(25))

# RNA Feature vs % MT UMI
FeatureScatter(gem5.raw, feature1 = "nFeature_RNA", feature2 = "percent.mt", log = FALSE) & 
  geom_hline(yintercept = c(0,25)) & 
  geom_vline(xintercept = c(450,4000))

# RNA Feature vs % MT UMI
FeatureScatter(gem5.raw, feature1 = "nCount_RNA", feature2 = "percent.mt", log = FALSE) & 
  geom_hline(yintercept = c(0,25)) & 
  geom_vline(xintercept = c(600,12500))

# 5) Calculate ATAC QC metrics ----

DefaultAssay(gem5.raw) <- "ATAC"

summary(gem5.raw$nCount_ATAC) # Median peaks - 
summary(gem5.raw$nFeature_ATAC) # Median features - 

# Calculation fraction of ATAC reads in peaks (FRiP)
gem5.raw$frac_reads_in_peaks <- gem5.raw$atac_peak_region_fragments / gem5.raw$atac_fragments
summary(gem5.raw$frac_reads_in_peaks) # FRiP median - 

# Calculate ATAC genomic blacklist fraction
gem5.raw$blacklist_fraction <- FractionCountsInRegion(
  object = gem5.raw, 
  assay = 'ATAC',
  regions = blacklist_hg38_unified)
summary(gem5.raw$blacklist_fraction) # Median blacklist fraction -

# Calculate transcription site enrichment score
gem5.raw <- TSSEnrichment(object = gem5.raw, fast = FALSE)
summary(gem5.raw$TSS.enrichment) # TSS enrichment median - 
gem5.raw$high.tss <- ifelse(gem5.raw$TSS.enrichment > 4, 'High', 'Low') # Split cells into groups based on chosen QC threshold

# Calculate nucleosome signal
gem5.raw <- NucleosomeSignal(object = gem5.raw)
summary(gem5.raw$nucleosome_signal) # Median nucleosome signal - 
gem5.raw$nucleosome_group <- ifelse(gem5.raw$nucleosome_signal > 0.7, 'NS > 0.7', 'NS < 0.7') # Split cells into groups based on chosen QC threshold

# 6) Select ATAC QC thresholds ----

# ATAC UMI Count
VlnPlot(gem5.raw, features = "nCount_ATAC", pt.size = 0.1, log = TRUE) + NoLegend() & geom_hline(yintercept = c(700,40000))
VlnPlot(gem5.raw, features = "nCount_ATAC", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(700,40000))
RidgePlot(gem5.raw, features = "nCount_ATAC", log = TRUE) + NoLegend() & geom_vline(xintercept = c(700,40000))
RidgePlot(gem5.raw, features = "nCount_ATAC", log = FALSE) + NoLegend() & geom_vline(xintercept = c(700,40000))

# ATAC Feature Count
VlnPlot(gem5.raw, features = "nFeature_ATAC", pt.size = 0.1, log = TRUE, , y.max = 15100) + NoLegend() & geom_hline(yintercept = c(350,15000))
VlnPlot(gem5.raw, features = "nFeature_ATAC", pt.size = 0.1, log = FALSE, y.max = 15100) + NoLegend() & geom_hline(yintercept = c(350,15000))
RidgePlot(gem5.raw, features = "nFeature_ATAC", log = TRUE) + NoLegend() & geom_vline(xintercept = c(350,15000))
RidgePlot(gem5.raw, features = "nFeature_ATAC", log = FALSE) + NoLegend() & geom_vline(xintercept = c(350,15000))

# ATAC Feature vs Count
FeatureScatter(gem5.raw, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC", log = TRUE) & 
  geom_hline(yintercept = c(350,15000)) & 
  geom_vline(xintercept = c(700,40000))
DensityScatter(gem5.raw, x = "nCount_ATAC", y = "nFeature_ATAC", log_x = TRUE, log_y = TRUE) & 
  geom_hline(yintercept = c(350,15000)) & 
  geom_vline(xintercept = c(700,40000))

# ATAC Peaks vs % MT UMI
FeatureScatter(gem5.raw, feature1 = "nCount_ATAC", feature2 = "percent.mt", log = FALSE) & 
  geom_hline(yintercept = c(0,25)) & 
  geom_vline(xintercept = c(700,40000))

# ATAC Features vs % MT UMI
FeatureScatter(gem5.raw, feature1 = "nFeature_ATAC", feature2 = "percent.mt", log = FALSE) & 
  geom_hline(yintercept = c(0,25)) & 
  geom_vline(xintercept = c(350,15000))

# ATAC Mitochondrial reads
VlnPlot(gem5.raw, features = "atac_mitochondrial_reads", pt.size = 0.1, log = TRUE) + NoLegend() & geom_hline(yintercept = c(0,5000))
VlnPlot(gem5.raw, features = "atac_mitochondrial_reads", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(0,5000))

# ATAC TSS Enrichment Score
VlnPlot(gem5.raw, features = "TSS.enrichment", pt.size = 0.1, log = TRUE, y.max = 30) + NoLegend() & geom_hline(yintercept = c(4))
VlnPlot(gem5.raw, features = "TSS.enrichment", pt.size = 0.1, log = FALSE, y.max = 30) + NoLegend() & geom_hline(yintercept = c(4))
RidgePlot(gem5.raw, features = "TSS.enrichment", log = TRUE) + NoLegend() & geom_vline(xintercept = c(4))
RidgePlot(gem5.raw, features = "TSS.enrichment", log = FALSE) + NoLegend() & geom_vline(xintercept = c(4))

# TSS enrichment plot grouped by high versus low enrichment score
TSSPlot(gem5.raw, group.by = 'high.tss')

# ATAC Nucleosome Signal
VlnPlot(gem5.raw, features = "nucleosome_signal", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(0.7))
RidgePlot(gem5.raw, features = "nucleosome_signal") + NoLegend() & geom_vline(xintercept = c(0.7))

# ATAC Nucleosome Signal vs TSS Enrichment Score
FeatureScatter(gem5.raw, feature1 = "TSS.enrichment", feature2 = "nucleosome_signal", log = TRUE) & 
  geom_hline(yintercept = c(0.7)) & 
  geom_vline(xintercept = c(4))
DensityScatter(gem5.raw, x = "TSS.enrichment", y = "nucleosome_signal", log_x = TRUE, log_y = TRUE) & 
  geom_hline(yintercept = c(0.75)) & 
  geom_vline(xintercept = c(4))

# Fragment length histogram grouped by high versus low nucleosome signal
FragmentHistogram(object = gem5.raw, group.by = 'nucleosome_group')

# ATAC FRiP
VlnPlot(gem5.raw, features = "frac_reads_in_peaks", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(0.5))
RidgePlot(gem5.raw, features = "frac_reads_in_peaks") + NoLegend() & geom_vline(xintercept = c(0.5))

# ATAC Blacklist Fraction
VlnPlot(gem5.raw, features = "blacklist_fraction", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(0.005))
RidgePlot(gem5.raw, features = "blacklist_fraction") + NoLegend() & geom_vline(xintercept = c(0.005))

# 7) Calculate ADT QC metrics ----

summary(gem5.raw$nCount_ADT)
summary(gem5.raw$nFeature_ADT)

# 8) Select ADT QC thresholds ----

# ADT Count
VlnPlot(gem5.raw, features = "nCount_ADT", pt.size = 0.1, log = TRUE) + NoLegend() & geom_hline(yintercept = c(175,2500))
VlnPlot(gem5.raw, features = "nCount_ADT", pt.size = 0.1, log = FALSE) + NoLegend() & geom_hline(yintercept = c(175,2500))
RidgePlot(gem5.raw, features = "nCount_ADT", log = TRUE) + NoLegend() & geom_vline(xintercept = c(175,2500))
RidgePlot(gem5.raw, features = "nCount_ADT", log = FALSE) + NoLegend() & geom_vline(xintercept = c(175,2500))

# 9) Calculate HTO QC metrics ----

summary(gem5.raw$nCount_HTO)

# 10) Select HTO QC thresholds ----

# HTO Count
VlnPlot(gem5.raw, features = "nCount_HTO", pt.size = 0.1, log = TRUE, y.max = 4100) + NoLegend() & geom_hline(yintercept = c(50,3000))
VlnPlot(gem5.raw, features = "nCount_HTO", pt.size = 0.1, log = FALSE, y.max = 4100) + NoLegend() & geom_hline(yintercept = c(50,3000))
RidgePlot(gem5.raw, features = "nCount_HTO", log = TRUE) + NoLegend() & geom_vline(xintercept = c(50,3000))

# 11) Save Seurat object with all QC metrics appended before filtering ----

gem5.prefilt <- gem5.raw
saveRDS(gem5.prefilt,"/filepath/step10_quality_control/gem5.prefilt.rds")

# 12) Filter Seurat object based on QC thresholds selected above for each modality ----

gem5.qc <- subset(x = gem5.prefilt, subset = 
                             nCount_RNA > 600 &
                             nCount_RNA < 12500 &
                             nFeature_RNA > 450 &
                             nFeature_RNA < 4000 &
                             percent.mt < 25 &
                             nCount_ATAC < 40000 &
                             nCount_ATAC > 700 &
                             nFeature_ATAC < 15000 &
                             nFeature_ATAC > 350 &
                             TSS.enrichment > 4 &
                             nucleosome_signal < 0.7 &
                             blacklist_fraction < 0.005 &
                             frac_reads_in_peaks > 0.5 &
                             atac_mitochondrial_reads < 5000 &
                             nCount_ADT > 175 &
                             nCount_ADT < 2500 &
                             nCount_HTO > 50 &
                             nCount_HTO < 3000 )

# 13) Assess number of cells filtered per sample ----

(ncol(gem5.prefilt))
(ncol(gem5.qc) - ncol(gem5.prefilt))
((ncol(gem5.qc) - ncol(gem5.prefilt)) / ncol(gem5.prefilt) * 100)
(ncol(gem5.qc))

# Extract metadata from Seurat objects, pre vs post QC
metadata_preQC <- gem5.prefilt@meta.data
metadata_postQC <- gem5.qc@meta.data

# Summarize the number of cells per sample, with QC "stage" set as "Before" vs "Filtered"
cell_counts_preQC <- metadata_preQC %>% 
  group_by(hto.sort) %>% 
  summarize(cell_count = n()) %>% 
  mutate(stage = 'Before QC')
cell_counts_postQC <- metadata_postQC %>% 
  group_by(hto.sort) %>% 
  summarize(cell_count = n()) %>% 
  mutate(stage = 'Filtered')

# Calculate total cell count, pre vs post QC
total_preQC <- data.frame(hto.sort = "All Samples", cell_count = sum(cell_counts_preQC$cell_count), stage = "Before QC")
total_postQC <- data.frame(hto.sort = "All Samples", cell_count = sum(cell_counts_postQC$cell_count), stage = "Filtered")

# Combine individual donor and overall cell counts into a single data frame
cell_counts_combined <- bind_rows(cell_counts_preQC, cell_counts_postQC, total_preQC, total_postQC)

# Calculate overall frequency
total_cells_preQC <- sum(cell_counts_preQC$cell_count)
total_cells_postQC <- sum(cell_counts_postQC$cell_count)

cell_counts_combined <- cell_counts_combined %>%
  mutate(frequency = ifelse(stage == 'Before QC', 
                            cell_count / total_cells_preQC, 
                            cell_count / total_cells_postQC))

# Define custom colors for the stages
colors <- c('Before QC' = 'lightgrey', 'Filtered' = 'aquamarine4')

# Create a bar plot with annotations
ggplot(cell_counts_combined, aes(x = hto.sort, y = cell_count, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = paste0(cell_count, "\n", round(frequency * 100, 1), "%")), 
            position = position_dodge(width = 0.97), vjust = -0.5,
            size = 3) +
  scale_fill_manual(values = colors) +
  labs(title = "Cell Number & Percent of Total - Pre vs Post QC (GEM5)",
       x = "Sample",
       y = "Number of Cells") +
  theme_minimal() +
  ylim(0,5000)

# 14) Save new object after QC filtering ----
saveRDS(gem5.qc, "/filepath/step10_quality_control/gem5.qc.rds")

# 15) Repeat for all GEM wells and proceed to next step ----