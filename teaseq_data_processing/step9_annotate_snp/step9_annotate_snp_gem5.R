# Step 9 - Annotate SNP Profile (GEM5) ----

# Matching unbiased each unbiased SNP cluster to expected HTO pattern per donor

# Substeps
# 1) Set up R working environment
# 2) Import singlet Seurat object for chosen GEM well
# 3) Visualize signal for each HTO across SNP profiles
# 4) Label donor of origin based on HTO signal for each SNP profile
# 5) Visualize both expected HTO signals for each annotated donor
# 6) Save donor-labeled Seurat object
# 7) Repeat for all GEM wells then proceed to next step

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
setwd("/filepath/step9_annotate_snp")

# 2) Import singlets Seurat object for chosen GEM well ----
gem5.singlet <- readRDS("/filepath/step8_filter_multiplets/gem5.singlet.umi.rds")

# 3) Visualize signal for each unique HTO across SNP profiles ----

# Set idents to SNP-based souporcell profiles
Idents(gem5.singlet) <- "snp.profile"

# CLR-normalize HTO counts across cells
DefaultAssay(gem5.singlet) <- "HTO"
gem5.singlet <- NormalizeData(gem5.singlet, assay= "HTO", normalization.method = "CLR", margin = 2)

# Plot HTO counts for each SNP profile
VlnPlot(gem5.singlet, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "snp.profile", pt.size = 0.1) + NoLegend()

# 4) Label donor of origin based on distinct HTO signal in each SNP profile ----

# Create cell identifier for donor based on HTO signal in each SNP profile
gem5.singlet$snp.donor <- gem5.singlet$snp.profile
Idents(gem5.singlet) <- "snp.donor"
gem5.annot <- RenameIdents(gem5.singlet,
                           "5" = "P1",
                           "4" = "P2",
                           "2" = "P3",
                           "0" = "P4",
                           "6" = "T1",
                           "3" = "T2",
                           "1" = "T3",
                           "7" = "T4")
gem5.annot$snp.donor <- Idents(gem5.annot)

# Plot HTO counts for each annotated SNP donor
VlnPlot(gem5.annot, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "snp.donor", pt.size = 0.1) + NoLegend()

# Verify HTO pattern between snp.donor and hto.donor groups match
VlnPlot(gem5.annot, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "hto.donor", pt.size = 0.1) + NoLegend()

# Create cell identifier for tissue based on HTO signal in each SNP profile
Idents(gem5.annot) <- "snp.donor"
gem5.annot$snp.tissue <- gem5.annot$snp.donor
Idents(gem5.annot) <- "snp.tissue"
gem5.annot <- RenameIdents(gem5.annot, 
                            "P1" = "PBMC",
                            "P2" = "PBMC",
                            "P3" = "PBMC",
                            "P4" = "PBMC",
                            "T1" = "Tonsil",
                            "T2" = "Tonsil",
                            "T3" = "Tonsil",
                            "T4" = "Tonsil")
gem5.annot$snp.tissue <- Idents(gem5.annot)

# Plot HTO counts for each annotated SNP tissue
VlnPlot(gem5.annot, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "snp.tissue", pt.size = 0.1) + NoLegend()

# Verify HTO pattern between hto.tissue and snp.tissue groups match
VlnPlot(gem5.annot, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "hto.tissue", pt.size = 0.1) + NoLegend()

# Check HTO expression pattern in HTO-defined flow sorting groups
VlnPlot(gem5.annot, features = c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10","HTO-11","HTO-12","HTO-13","HTO-14","HTO-15","HTO-16"), 
        ncol = 4, group.by = "hto.sort", pt.size = 0.1) + NoLegend()

# 5) Visualize both expected HTO signals for each annotated donor by SNP ----

# Visualize signal for HTO 1 & 5 in Donor P1
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P1")
ggplot(hto_df, aes(x = `HTO-1`, y = `HTO-5`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-1", y = "HTO-5") +
  theme_minimal() +
  ggtitle("HTO-1 and HTO-5 Signal in Donor P1") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 2 & 6 in Donor P2
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P2")
ggplot(hto_df, aes(x = `HTO-2`, y = `HTO-6`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-2", y = "HTO-6") +
  theme_minimal() +
  ggtitle("HTO-2 and HTO-6 Signal in Donor P2") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 3 & 7 in Donor P3
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P3")
ggplot(hto_df, aes(x = `HTO-3`, y = `HTO-7`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-3", y = "HTO-7") +
  theme_minimal() +
  ggtitle("HTO-3 and HTO-7 Signal in Donor P3") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 4 & 8 in Donor P4
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "P4")
ggplot(hto_df, aes(x = `HTO-4`, y = `HTO-8`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-4", y = "HTO-8") +
  theme_minimal() +
  ggtitle("HTO-4 and HTO-8 Signal in Donor P4") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 9 & 13 in Donor T1
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T1")
ggplot(hto_df, aes(x = `HTO-9`, y = `HTO-13`)) +
  geom_point() +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-9", y = "HTO-13") +
  theme_minimal() +
  ggtitle("HTO-9 and HTO-13 Signal in Donor T1") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 10 & 14 in Donor T2
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T2")
ggplot(hto_df, aes(x = `HTO-10`, y = `HTO-14`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-10", y = "HTO-14") +
  theme_minimal() +
  ggtitle("HTO-10 and HTO-14 Signal in Donor T2") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 11 & 15 in Donor T3
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T3")
ggplot(hto_df, aes(x = `HTO-11`, y = `HTO-15`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-11", y = "HTO-15") +
  theme_minimal() +
  ggtitle("HTO-11 and HTO-15 Signal in Donor T3") +
  guides(fill=guide_legend(title="ndensity"))

# Visualize signal for HTO 12 & 16 in Donor T4
hto_df <- t(as.data.frame(GetAssayData(subset(gem5.annot), assay = "HTO")))
snp <- as.data.frame(gem5.annot$snp.donor)
hto_df <- cbind(hto_df, snp)
colnames(hto_df)[17] <- "SNP"
hto_df <- subset(hto_df, SNP == "T4")
ggplot(hto_df, aes(x = `HTO-12`, y = `HTO-16`)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 0.25, colour = "black") +
  geom_density_2d_filled(alpha = 0.5, contour_var = "ndensity") +
  labs(x = "HTO-12", y = "HTO-16") +
  theme_minimal() +
  ggtitle("HTO-12 and HTO-16 Signal in Donor T4") +
  guides(fill=guide_legend(title="ndensity"))

# 6) Save SNP-annotated Seurat object ----
saveRDS(gem5.annot, "/filepath/step9_annotate_snp/gem5.annot.snp.rds")

# 7) Repeat for all GEM wells then proceed to next step ----