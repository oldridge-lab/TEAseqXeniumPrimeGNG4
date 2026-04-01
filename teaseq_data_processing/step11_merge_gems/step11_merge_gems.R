# Step 11 - Create unified ATAC peak matrix across all eight GEM wells and merge Seurat objects  ----

# Substeps
# 1) Set up R working environment
# 2) Import QC'd Seurat objects for each GEM well from Step 10
# 3) Unify ATAC peaks between GEM wells that were previously called by cellranger-arc independently
# 4) Create and assign unified ATAC peak count matrix for each GEM well
# 5) Merge Seurat objects from each GEM well
# 6) Join RNA assay count layers between original GEM wells 
# 7) Organize merged object annotations
# 8) Verify expected HTO pattern in annotated sorting groups, donors, and tissues
# 9) Assess number of cells filtered per sample
# 10) Save merged Seurat object
# 11) Proceed to next step - cell-cycle gene scoring

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
library(viridis) # 0.6.5

# Set working directory for this step
setwd("/filepath/step11_merge_gems")

# 2) Import Seurat objects for each GEM well (empty-filtered, multiplet-filtered, HTO-demultiplexed, and SNP-demultiplexed) ----
gem1.obj <- readRDS("/filepath/step10_quality_control/gem1.qc.rds") 
gem2.obj <- readRDS("/filepath/step10_quality_control/gem2.qc.rds") 
gem3.obj <- readRDS("/filepath/step10_quality_control/gem3.qc.rds") 
gem4.obj <- readRDS("/filepath/step10_quality_control/gem4.qc.rds") 
gem5.obj <- readRDS("/filepath/step10_quality_control/gem5.qc.rds") 
gem6.obj <- readRDS("/filepath/step10_quality_control/gem6.qc.rds") 
gem7.obj <- readRDS("/filepath/step10_quality_control/gem7.qc.rds") 
gem8.obj <- readRDS("/filepath/step10_quality_control/gem8.qc.rds") 

# Check cell number for each GEM well
(gem1.n <- ncol(gem1.obj)) # 3806
(gem2.n <- ncol(gem2.obj)) # 3949
(gem3.n <- ncol(gem3.obj)) # 4303
(gem4.n <- ncol(gem4.obj)) # 3732
(gem5.n <- ncol(gem5.obj)) # 4381
(gem6.n <- ncol(gem6.obj)) # 4263
(gem7.n <- ncol(gem7.obj)) # 3942
(gem8.n <- ncol(gem8.obj)) # 3830

# Expected cell number after merging GEM wells
(merged.n <- gem1.n + gem2.n + gem3.n + gem4.n + gem5.n + gem6.n + gem7.n + gem8.n) # 32206 cells

# 3) Unify ATAC peaks between GEM wells that were previously called by cellranger-arc independently ----
combined.peaks <- UnifyPeaks(object.list = list(gem1.obj[['ATAC']], 
                                                gem2.obj[['ATAC']],
                                                gem3.obj[['ATAC']],
                                                gem4.obj[['ATAC']],
                                                gem5.obj[['ATAC']],
                                                gem6.obj[['ATAC']],
                                                gem7.obj[['ATAC']],
                                                gem8.obj[['ATAC']])
                             , mode = "reduce")
peakwidths <- width(combined.peaks)
length(peakwidths)
summary(peakwidths)

# Save combined peak file
saveRDS(combined.peaks, '/filepath/step11_merge_gems/combined.peaks.rds')

# 4) Create and assign unified ATAC peak count matrix for each GEM well ----

# Retrieve Genome Reference Consortium Human Build 38 (GRCh38) annotation file previously generated during Step 5 - Preprocess GEMs
annotations <- readRDS('/filepath/step5_preprocess_gems/atac_annotations.rds')

# GEM 1
gem1_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem1.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem1.obj)
)
frag_1.file <- "/filepath/cellranger-arc_output/GEM1_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem1.obj) <- 'RNA'
gem1.obj[['ATAC']] <- NULL 
gem1.obj[['ATAC']] <- CreateChromatinAssay(counts = gem1_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_1.file,
                                               annotation = annotations)
saveRDS(gem1_unified.counts, '/filepath/step11_merge_gems/gem1_unified.counts.rds')

# GEM 2
gem2_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem2.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem2.obj)
)
frag_2.file <- "/filepath/cellranger-arc_output/GEM2_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem2.obj) <- 'RNA'
gem2.obj[['ATAC']] <- NULL
gem2.obj[['ATAC']] <- CreateChromatinAssay(counts = gem2_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_2.file,
                                               annotation = annotations)
saveRDS(gem2_unified.counts, '/filepath/step11_merge_gems/gem2_unified.counts.rds')

# GEM 3
gem3_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem3.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem3.obj)
)
frag_3.file <- "/filepath/cellranger-arc_output/GEM3_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem3.obj) <- 'RNA'
gem3.obj[['ATAC']] <- NULL
gem3.obj[['ATAC']] <- CreateChromatinAssay(counts = gem3_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_3.file,
                                               annotation = annotations)
saveRDS(gem3_unified.counts, '/filepath/step11_merge_gems/gem3_unified.counts.rds')

# GEM 4
gem4_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem4.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem4.obj)
)
frag_4.file <- "/filepath/cellranger-arc_output/GEM4_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem4.obj) <- 'RNA'
gem4.obj[['ATAC']] <- NULL
gem4.obj[['ATAC']] <- CreateChromatinAssay(counts = gem4_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_4.file,
                                               annotation = annotations)
saveRDS(gem4_unified.counts, '/filepath/step11_merge_gems/gem4_unified.counts.rds')

# GEM 5
gem5_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem5.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem5.obj)
)
frag_5.file <- "/filepath/cellranger-arc_output/GEM5_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem5.obj) <- 'RNA'
gem5.obj[['ATAC']] <- NULL
gem5.obj[['ATAC']] <- CreateChromatinAssay(counts = gem5_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_5.file,
                                               annotation = annotations)
saveRDS(gem5_unified.counts, '/filepath/step11_merge_gems/gem5_unified.counts.rds')

# GEM 6
gem6_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem6.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem6.obj)
)
frag_6.file <- "/filepath/cellranger-arc_output/GEM6_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem6.obj) <- 'RNA'
gem6.obj[['ATAC']] <- NULL
gem6.obj[['ATAC']] <- CreateChromatinAssay(counts = gem6_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_6.file,
                                               annotation = annotations)
saveRDS(gem6_unified.counts, '/filepath/step11_merge_gems/gem6_unified.counts.rds')

# GEM 7
gem7_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem7.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem7.obj)
)
frag_7.file <- "/filepath/cellranger-arc_output/GEM7_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem7.obj) <- 'RNA'
gem7.obj[['ATAC']] <- NULL
gem7.obj[['ATAC']] <- CreateChromatinAssay(counts = gem7_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_7.file,
                                               annotation = annotations)
saveRDS(gem7_unified.counts, '/filepath/step11_merge_gems/gem7_unified.counts.rds')

# GEM 8
gem8_unified.counts <- FeatureMatrix(
  fragments = Fragments(gem8.obj[['ATAC']]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(gem8.obj)
)
frag_8.file <- "/filepath/cellranger-arc_output/GEM8_ARC/outs/atac_fragments.tsv.gz"
DefaultAssay(gem8.obj) <- 'RNA'
gem8.obj[['ATAC']] <- NULL
gem8.obj[['ATAC']] <- CreateChromatinAssay(counts = gem8_unified.counts,
                                               sep = c(":", "-"),
                                               fragments = frag_8.file,
                                               annotation = annotations)
saveRDS(gem8_unified.counts, '/filepath/step11_merge_gems/gem8_unified.counts.rds')

# 5) Merge GEM well Seurat objects ----
merged.obj <- merge(gem1.obj, 
                    y = c(gem2.obj, 
                          gem3.obj,
                          gem4.obj,
                          gem5.obj,
                          gem6.obj,
                          gem7.obj,
                          gem8.obj),
                    add.cell.ids = c("GEM1",
                                     "GEM2",
                                     "GEM3",
                                     "GEM4",
                                     "GEM5",
                                     "GEM6",
                                     "GEM7",
                                     "GEM8"), 
                    project = "GEM_merged")

# Check that merged object barcode number is equivalent to sum of independent GEM objects
ncol(merged.obj) == merged.n

# 6) Join RNA assay count layers between original GEM wells ----

DefaultAssay(merged.obj) <- "RNA"
merged.obj <- JoinLayers(merged.obj)

# 7) Organize merged object annotations ----

merged.obj$hto.donor <- factor(x = merged.obj$hto.donor, levels = c('P1','P2','P3','P4','T1','T2','T3','T4'))
merged.obj$snp.donor <- factor(x = merged.obj$snp.donor, levels = c('P1','P2','P3','P4','T1','T2','T3','T4'))
merged.obj$hto.max <- factor(x = merged.obj$hto.max, levels = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"))
merged.obj$hto.sort <- factor(x = merged.obj$hto.sort, levels = c("P1-Bulk","P2-Bulk","P3-Bulk","P4-Bulk","P1-CD4","P2-CD4","P3-CD4","P4-CD4","T1-Bulk","T2-Bulk","T3-Bulk","T4-Bulk","T1-CD4","T2-CD4","T3-CD4","T4-CD4"))

# Create 'assigned sex at birth' (asab) label for object
merged.obj$asab <- merged.obj$hto.donor
Idents(merged.obj) <- "asab"
merged.obj <- RenameIdents(merged.obj,
                                 "P1" = "Female",
                                 "P2" = "Female",
                                 "P3" = "Male",
                                 "P4" = "Male",
                                 "T1" = "Male",
                                 "T2" = "Male",
                                 "T3" = "Female",
                                 "T4" = "Female")
merged.obj$asab <- Idents(merged.obj)
Idents(merged.obj) <- "hto.donor"

# 8) Verify expected HTO pattern in annotated sorting groups, donors, and tissues ----

DefaultAssay(merged.obj) <- 'HTO'
merged.obj <- NormalizeData(merged.obj, margin = 2, normalization.method = 'CLR', assay = 'HTO')

# Plot HTO counts for each annotated SNP donor
VlnPlot(merged.obj, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "snp.donor", pt.size = 0.1) + NoLegend()

# Verify HTO pattern between snp.donor and hto.donor groups match
VlnPlot(merged.obj, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "hto.donor", pt.size = 0.1) + NoLegend()

# Plot HTO counts for each annotated SNP tissue
VlnPlot(merged.obj, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "snp.tissue", pt.size = 0.1) + NoLegend()

# Verify HTO pattern between hto.tissue and snp.tissue groups match
VlnPlot(merged.obj, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "hto.tissue", pt.size = 0.1) + NoLegend()

# Check HTO expression pattern in HTO-defined flow sorting groups
VlnPlot(merged.obj, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "hto.sort", pt.size = 0.1) + NoLegend()

# Visualize signal for HTO 1 & 5 in Donor P1
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P1")
ggplot(hto_df, aes(x = `HTO-1`, y = `HTO-5`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-1", y = "HTO-5") +
  theme_minimal() +
  ggtitle("HTO-1 and HTO-5 Signal in Donor P1") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 2 & 6 in Donor P2
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P2")
ggplot(hto_df, aes(x = `HTO-2`, y = `HTO-6`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-2", y = "HTO-6") +
  theme_minimal() +
  ggtitle("HTO-2 and HTO-6 Signal in Donor P2") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 3 & 7 in Donor P3
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P3")
ggplot(hto_df, aes(x = `HTO-3`, y = `HTO-7`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-3", y = "HTO-7") +
  theme_minimal() +
  ggtitle("HTO-3 and HTO-7 Signal in Donor P3") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 4 & 8 in Donor P4
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P4")
ggplot(hto_df, aes(x = `HTO-4`, y = `HTO-8`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-4", y = "HTO-8") +
  theme_minimal() +
  ggtitle("HTO-4 and HTO-8 Signal in Donor P4") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 9 & 13 in Donor T1
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T1")
ggplot(hto_df, aes(x = `HTO-9`, y = `HTO-13`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-9", y = "HTO-13") +
  theme_minimal() +
  ggtitle("HTO-9 and HTO-13 Signal in Donor T1") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 10 & 14 in Donor T2
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T2")
ggplot(hto_df, aes(x = `HTO-10`, y = `HTO-14`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-10", y = "HTO-14") +
  theme_minimal() +
  ggtitle("HTO-10 and HTO-14 Signal in Donor T2") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 11 & 15 in Donor T3
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T3")
ggplot(hto_df, aes(x = `HTO-11`, y = `HTO-15`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-11", y = "HTO-15") +
  theme_minimal() +
  ggtitle("HTO-11 and HTO-15 Signal in Donor T3") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 12 & 16 in Donor T4
hto_df <- t(as.data.frame(GetAssayData(subset(merged.obj), assay = "HTO")))
snp <- as.data.frame(merged.obj$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T4")
ggplot(hto_df, aes(x = `HTO-12`, y = `HTO-16`)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-12", y = "HTO-16") +
  theme_minimal() +
  ggtitle("HTO-12 and HTO-16 Signal in Donor T4") +
  guides(fill=guide_legend(title="ndensity"))

# Compare tissue annotations based on HTO and SNP pattern
annot.df <- data.frame(merged.obj$snp.tissue, merged.obj$hto.tissue)
colnames(annot.df) <- c('SNP','HTO')
levels <- c('PBMC','Tonsil')
annot.df$SNP <- factor(annot.df$SNP, levels = levels)
annot.df$HTO <- factor(annot.df$HTO, levels = levels)
conf_matrix <- confusionMatrix(annot.df$SNP, annot.df$HTO)
conf_table <- as.data.frame(conf_matrix$table)
colnames(conf_table) <- c("SNP", "HTO", "Cell Number")
ggplot(conf_table, aes(x = SNP, y = HTO, fill = `Cell Number`)) +
  geom_tile() +
  geom_text(aes(label = `Cell Number`), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "coral") +
  theme_minimal() +
  labs(title = "HTO vs SNP Tissue Annotation",
       x = "SNP Assignment",
       y = "HTO Assignment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compare donor annotations based on HTO and SNP pattern
annot.df <- data.frame(merged.obj$snp.donor, merged.obj$hto.donor)
colnames(annot.df) <- c('SNP','HTO')
conf_matrix <- confusionMatrix(annot.df$SNP, annot.df$HTO)
conf_table <- as.data.frame(conf_matrix$table)
colnames(conf_table) <- c("SNP", "HTO", "Cell Number")
ggplot(conf_table, aes(x = SNP, y = HTO, fill = `Cell Number`)) +
  geom_tile() +
  geom_text(aes(label = `Cell Number`), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "coral") +
  theme_minimal() +
  labs(title = "HTO vs SNP Donor Annotation",
       x = "SNP Assignment",
       y = "HTO Assignment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 9) Assess number of cells filtered per sample ----

table(merged.obj$hto.sort)

# Extract metadata from object
merged.md <- merged.obj@meta.data

# Summarize the number of cells per sample
cell_counts <- merged.md %>% group_by(hto.sort) %>% summarize(cell_count = n())

# Calculate total cell count, pre vs post QC
total_count <- sum(cell_counts$cell_count)

# Calculate overall frequency
cell_counts <- cell_counts %>% mutate(frequency = cell_count / total_count)

# Remove All Samples column
cell_counts <- cell_counts[1:16,]

# Create a bar plot
ggplot(cell_counts, aes(x = hto.sort, y = cell_count, fill = frequency)) + geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = paste0(cell_count, "\n", round(frequency * 100, 2), "%")), vjust = -0.5, size = 5) +
  labs(title = "FACS Sample Composition - Cell Number and Frequency of Pool Total after QC",
       x = "FACS Sample",
       y = "Cell Number") +
  theme_minimal() +
  scale_fill_viridis() +
  ylim(0,3750) +
  theme(
    plot.title = element_text(size = 25),    
    axis.title.x = element_text(size = 25),   
    axis.title.y = element_text(size = 25, margin = margin(r = 15)),  
    axis.text.x = element_text(size = 15, margin = margin(b = 15)),   
    axis.text.y = element_text(size = 10)     
  )

# 10) Save merged Seurat object ----

saveRDS(merged.obj, "/filepath/step11_merge_gems/merged.obj.rds")

# 11) Proceed to next step - cell-cycle gene scoring ----