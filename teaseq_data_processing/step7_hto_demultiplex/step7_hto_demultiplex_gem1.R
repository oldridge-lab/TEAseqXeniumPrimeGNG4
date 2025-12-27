# Step 7 - HTO Demultiplex (GEM1)

# Calculating % HTO UMI metrics

# Substeps ----

# 1) Set up R working environment
# 2) Import Seurat object with empty droplet removed for GEM well
# 3) Calculate % HTO UMI metrics
# 4) Annotate donor, sort group, and tissue based on % HTO UMI metrics
# 5) Save HTO-demultiplexed Seurat object
# 6) Repeat for all 8 GEM wells and proceed to next step

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
setwd("/filepath/step7_hto_demultiplex")

# 2) Import Seurat object with empty droplets removed for GEM well ----
gem1.filt <- readRDS('/filepath/step6_remove_empty/gem1.filt.rds')

# 3) Calculate % HTO UMI metrics ----

# Calculate % HTO UMI metrics
hto.counts.mtx <- gem1.filt[['HTO']]$counts
hto.counts.mtx.sort <- apply(hto.counts.mtx, 2, function(x) sort(x, decreasing = TRUE))
hto.counts.mtx.sort <- t(hto.counts.mtx.sort)
hto.counts.mtx.total <- apply(hto.counts.mtx.sort, 1, sum)
hto.counts.mtx.sort.pct <- hto.counts.mtx.sort / hto.counts.mtx.total * 100

# Append % HTO UMI metrics to Seurat object
gem1.filt$pct.hto.max <- hto.counts.mtx.sort.pct[,1]
gem1.filt$pct.hto.2nd <- hto.counts.mtx.sort.pct[,2]
gem1.filt$pct.hto.bg <- rowSums(hto.counts.mtx.sort.pct[,2:16])

# Find HTO with maximum % UMI in each cell 
hto.max.index <- apply(hto.counts.mtx, 2, function(col) {which.max(col)})
hto.max <- paste("HTO-", hto.max.index, sep = "")
gem1.filt$hto.max <- hto.max

# 5) Annotate cells based on % HTO UMI metrics ----

# Assign flow sorting group labels from % Max HTO metrics
gem1.filt$hto.sort <- gem1.filt$hto.max
Idents(gem1.filt) <- 'hto.sort'
gem1.filt <- RenameIdents(gem1.filt,
              "HTO-1" = "P1-Bulk", 
              "HTO-2" = "P2-Bulk", 
              "HTO-3" = "P3-Bulk", 
              "HTO-4" = "P4-Bulk", 
              "HTO-5" = "P1-CD4", 
              "HTO-6" = "P2-CD4",
              "HTO-7" = "P3-CD4",
              "HTO-8" = "P4-CD4", 
              "HTO-9" = "T1-Bulk",
              "HTO-10" = "T2-Bulk", 
              "HTO-11" = "T3-Bulk",
              "HTO-12" = "T4-Bulk",
              "HTO-13" = "T1-CD4",
              "HTO-14" = "T2-CD4",
              "HTO-15" = "T3-CD4", 
              "HTO-16" = "T4-CD4")
gem1.filt$hto.sort <- Idents(gem1.filt)
as.data.frame(gem1.filt$hto.sort) %>% 
  ggplot(aes(x = gem1.filt$hto.sort)) + 
  geom_bar(fill = 'azure4') + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'Cells per flow sorting group (GEM1)',
       x = 'Flow sorting group',
       y = 'Number of cells') +
  theme_minimal()

# Assign donor labels from % Max HTO metrics
gem1.filt$hto.donor <- gem1.filt$hto.max
Idents(gem1.filt) <- 'hto.donor'
gem1.filt <- RenameIdents(gem1.filt,
                          "HTO-1" = "P1", 
                          "HTO-2" = "P2", 
                          "HTO-3" = "P3", 
                          "HTO-4" = "P4", 
                          "HTO-5" = "P1",
                          "HTO-6" = "P2", 
                          "HTO-7" = "P3",
                          "HTO-8" = "P4",
                          "HTO-9" = "T1", 
                          "HTO-10" = "T2",
                          "HTO-11" = "T3",
                          "HTO-12" = "T4", 
                          "HTO-13" = "T1",
                          "HTO-14" = "T2",
                          "HTO-15" = "T3", 
                          "HTO-16" = "T4")
gem1.filt$hto.donor <- Idents(gem1.filt)
as.data.frame(gem1.filt$hto.donor) %>% 
  ggplot(aes(x = gem1.filt$hto.donor)) + 
  geom_bar(fill = 'azure4') + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'Cells per donor (GEM1)',
       x = 'Donor',
       y = 'Number of cells') +
  theme_minimal()

# Assign tissue type labels from % Max HTO metrics
gem1.filt$hto.tissue <- gem1.filt$hto.max
Idents(gem1.filt) <- 'hto.tissue'
gem1.filt <- RenameIdents(gem1.filt,
                          "HTO-1" = "PBMC", 
                          "HTO-2" = "PBMC", 
                          "HTO-3" = "PBMC", 
                          "HTO-4" = "PBMC", 
                          "HTO-5" = "PBMC",
                          "HTO-6" = "PBMC", 
                          "HTO-7" = "PBMC",
                          "HTO-8" = "PBMC",
                          "HTO-9" = "Tonsil", 
                          "HTO-10" = "Tonsil",
                          "HTO-11" = "Tonsil",
                          "HTO-12" = "Tonsil", 
                          "HTO-13" = "Tonsil",
                          "HTO-14" = "Tonsil",
                          "HTO-15" = "Tonsil", 
                          "HTO-16" = "Tonsil")
gem1.filt$hto.tissue <- Idents(gem1.filt)
as.data.frame(gem1.filt$hto.tissue) %>% 
  ggplot(aes(x = gem1.filt$hto.tissue)) + 
  geom_bar(fill = 'azure4') + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'Cells per tissue (GEM1)',
       x = 'Tissue',
       y = 'Number of cells') +
  theme_minimal()

# Clean up excess or redundant metadata columns
gem1.filt$HTO_classification.global <- NULL
gem1.filt$HTO_classification <- NULL
gem1.filt$HTO_margin <- NULL
gem1.filt$HTO_maxID <- NULL
gem1.filt$hash.ID <- NULL
gem1.filt$HTO_secondID <- NULL
gem1.filt$orig.ident <- NULL

# 6) Save HTO-demultiplexed Seurat object ----
saveRDS(gem1.filt, '/filepath/step7_hto_demultiplex/gem1.hto.rds')

# 6) Repeat for all 8 GEM wells and proceed to next step ----