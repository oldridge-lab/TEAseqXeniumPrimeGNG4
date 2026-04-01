# Step 8 - Filter Multiplets (GEM1)

# Removing multiplets based on combined % HTO UMI, scDblFinder RNA/ATAC, and Souporcell methods

# Substeps ----

# 1) Set up R working environment
# 2) Import HTO-demultiplexed Seurat object for chosen GEM well
# 3) Run scDblFinder using RNA modality
# 4) Run scDblFinder using ATAC modality
# 5) Determine multiplet status from % HTO UMI metrics
# 6) Analyze scDblFinder RNA performance
# 7) Analyze scDblFinder ATAC performance
# 8) Analyze Souporcell performance
# 9) Analyze combined multiplet classifications
# 10) Filter multiplets identified using combined % HTO UMI, scDblFinder RNA/ATAC, and Souporcell approach - keep flags from other methods
# 11) Filter possibly missed multiplets with exceptionally high UMI counts
# 12) Save singlet Seurat object
# 13) Repeat for all 8 GEMs and proceed to next step

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
setwd("/filepath/step8_filter_multiplets")

# 2) Import HTO-demultiplexed Seurat object for GEM well ----
gem1.hto <- readRDS('/filepath/step7_hto_demultiplex/gem1.hto.rds')

# 3) Run scDblFinder using RNA modality ----

# Estimate expected multiplet rate using loaded cell number (20000) and 10x Multiome rate estimate table
# Chromium Next GEM Single Cell Multiome ATAC + Gene Expression Reagent Kits (CG000338 Rev F, August 26 2022)
(dbr.est <- 7.7/15400*20000*1/100)

# Convert Seurat object to SCE
DefaultAssay(gem1.hto) <- 'RNA'
gem1.sce <- as.SingleCellExperiment(gem1.hto)

# Run scDblFinder on RNA data
scDbl.rna_gem1 <- scDblFinder(gem1.sce, dbr = dbr.est) 

# Append scDblFinder outputs from RNA modality to original Seurat object
gem1.hto$scDblF.RNA.score <- scDbl.rna_gem1$scDblFinder.score
gem1.hto$scDblF.RNA.class <- scDbl.rna_gem1$scDblFinder.class
gem1.hto$scDblF.RNA.class <- ifelse(gem1.hto$scDblF.RNA.class == "singlet", "singlet", "multiplet")

# 4) Run scDblFinder using ATAC modality ----

DefaultAssay(gem1.hto) <- 'ATAC'
gem1.sce <- as.SingleCellExperiment(gem1.hto)
scDbl.atac_gem1 <- scDblFinder(gem1.sce, dbr = dbr.est, aggregateFeatures = TRUE, nfeatures = 25, processing = 'normFeatures') 

# Append scDblFinder outputs from ATAC modality to original Seurat object
gem1.hto$scDblF.ATAC.score <- scDbl.atac_gem1$scDblFinder.score
gem1.hto$scDblF.ATAC.class <- scDbl.atac_gem1$scDblFinder.class
gem1.hto$scDblF.ATAC.class <- ifelse(gem1.hto$scDblF.ATAC.class == "singlet", "singlet", "multiplet")

# 5) Determine multiplet status from % HTO UMI metrics ----

# Create dataframe containing % UMI from maximum, 2nd, and cumulative 2nd-16th HTO per cell
pct.hto.df <- data.frame(gem1.hto$pct.hto.max, gem1.hto$pct.hto.2nd, gem1.hto$pct.hto.bg)
colnames(pct.hto.df) <- c('pct.hto.max','pct.hto.2nd','pct.hto.bg')

# Assign multiplet status based on percentile thresholds
pct.hto.df <- pct.hto.df %>%
  mutate(hto.class = case_when(
    pct.hto.max > 50 & pct.hto.2nd < 15  & pct.hto.bg < 50 ~ 'singlet',
    pct.hto.2nd >=15 ~ 'multiplet',
    TRUE ~ 'negative'
  ))
pct.hto.df$hto.class <- factor(pct.hto.df$hto.class, levels = c('singlet','multiplet','negative'))
gem1.hto$hto.class <- pct.hto.df$hto.class

# Visualize results of % HTO UMI multiplet assignment
pct.hto.df %>% ggplot(aes(x = hto.class, fill = hto.class)) + 
  geom_bar() + 
  scale_fill_manual(values = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = '% HTO UMI Assignment (GEM1)',
       x = 'Assignment',
       y = 'Number of cells') +
  theme_minimal()

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by multiplet status
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = hto.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "hto.class")

# Visualize % UMI from maximum vs 2nd-16th HTO in each cell, colored by multiplet status
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.bg`, color = hto.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd-16th HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd-16th HTO (GEM1)") +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "Droplet")

# Compare UMI from each modality between % HTO UMI-based multiplet assignment classes
gem1.hto$hto.class <- pct.hto.df$hto.class
Idents(gem1.hto) <- 'hto.class'
gem1.hto@active.ident <- factor(gem1.hto@active.ident, levels = c('singlet','multiplet','negative'))
VlnPlot(gem1.hto, features = 'nCount_HTO', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ADT', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_RNA', pt.size = 0.1, log = TRUE,  cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()

# 6) Analyze scDblFinder RNA performance ----

# Add scDblFinder RNA classifications to % HTO UMI dataframe
pct.hto.df$scDblF.RNA.class <- gem1.hto$scDblF.RNA.class

# Visualize results of scDblF.RNA.class multiplet assignment
pct.hto.df %>% ggplot(aes(x = scDblF.RNA.class, fill = scDblF.RNA.class)) + 
  geom_bar() + 
  scale_fill_manual(values = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'scDblF.RNA.class Assignment (GEM1)',
       x = 'Assignment',
       y = 'Number of cells') +
  theme_minimal()

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by scDblF.RNA.class classification
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = scDblF.RNA.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO - colored by scDblF.RNA.class classification (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "scDblF.RNA.class")

# Visualize % UMI from maximum vs 2nd-16th HTO in each cell, colored by scDblF.RNA.class classification
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.bg`, color = scDblF.RNA.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd-16th HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd-16th HTO - colored by scDblF.RNA.class classification (GEM1)") +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "scDblF.RNA.class")

# Compare UMI from each modality between scDblF.RNA.class classes
Idents(gem1.hto) <- 'scDblF.RNA.class'
gem1.hto@active.ident <- factor(gem1.hto@active.ident, levels = c('singlet','multiplet','negative'))
VlnPlot(gem1.hto, features = 'nCount_HTO', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ADT', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_RNA', pt.size = 0.1, log = TRUE,  cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()

# 7) Analyze scDblFinder ATAC performance ----

# Add scDblFinder RNA classifications to % HTO UMI dataframe
pct.hto.df$scDblF.ATAC.class <- gem1.hto$scDblF.ATAC.class

# Visualize results of scDblF.ATAC.class multiplet assignment
pct.hto.df %>% ggplot(aes(x = scDblF.ATAC.class, fill = scDblF.ATAC.class)) + 
  geom_bar() + 
  scale_fill_manual(values = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'scDblF.ATAC.class Assignment (GEM1)',
       x = 'Assignment',
       y = 'Number of cells') +
  theme_minimal()

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by scDblF.ATAC.class classification
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = scDblF.ATAC.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO - colored by scDblF.ATAC.class classification (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "scDblF.ATAC.class")

# Visualize % UMI from maximum vs 2nd-16th HTO in each cell, colored by scDblF.ATAC.class classification
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.bg`, color = scDblF.ATAC.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd-16th HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd-16th HTO - colored by scDblF.ATAC.class classification (GEM1)") +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "scDblF.ATAC.class")

# Compare UMI from each modality between scDblF.ATAC.class classes
Idents(gem1.hto) <- 'scDblF.ATAC.class'
gem1.hto@active.ident <- factor(gem1.hto@active.ident, levels = c('singlet','multiplet','negative'))
VlnPlot(gem1.hto, features = 'nCount_HTO', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ADT', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_RNA', pt.size = 0.1, log = TRUE,  cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()

# 8) Analyze Souporcell performance ----

# Rename Souporcell metadata columns
gem1.hto$snp.class <- gem1.hto$status
gem1.hto@meta.data <- gem1.hto@meta.data %>%
  mutate(snp.class = recode(snp.class, 
                                 'unassigned' = 'negative', 
                                 'doublet' = 'multiplet'))
gem1.hto$snp.profile <- gem1.hto$assignment
gem1.hto$snp.singlet.logp <- gem1.hto$log_prob_singleton
gem1.hto$snp.multip.logp <- gem1.hto$log_prob_doublet
gem1.hto$snp.profile0.logp <- gem1.hto$cluster0
gem1.hto$snp.profile1.logp <- gem1.hto$cluster1
gem1.hto$snp.profile2.logp <- gem1.hto$cluster2
gem1.hto$snp.profile3.logp <- gem1.hto$cluster3
gem1.hto$snp.profile4.logp <- gem1.hto$cluster4
gem1.hto$snp.profile5.logp <- gem1.hto$cluster5
gem1.hto$snp.profile6.logp <- gem1.hto$cluster6
gem1.hto$snp.profile7.logp <- gem1.hto$cluster7

# Clean up excess Souporcell metadata
gem1.hto$status <- NULL
gem1.hto$assignment <- NULL
gem1.hto$log_prob_singleton <- NULL
gem1.hto$log_prob_doublet <- NULL
gem1.hto$cluster0 <- NULL
gem1.hto$cluster1 <- NULL
gem1.hto$cluster2 <- NULL
gem1.hto$cluster3 <- NULL
gem1.hto$cluster4 <- NULL
gem1.hto$cluster5 <- NULL
gem1.hto$cluster6 <- NULL
gem1.hto$cluster7 <- NULL

# Add Souporcell classifications to % HTO UMI dataframe
pct.hto.df$snp.class <- gem1.hto$snp.class

# Visualize results of Souporcell multiplet assignment
pct.hto.df %>% ggplot(aes(x = snp.class, fill = snp.class)) + 
  geom_bar() + 
  scale_fill_manual(values = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'Souporcell Assignment (GEM1)',
       x = 'Assignment',
       y = 'Number of cells') +
  theme_minimal()

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by Souporcell classification
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = snp.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO - colored by Souporcell classification (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "snp.class")

# Visualize % UMI from maximum vs 2nd-16th HTO in each cell, colored by Souporcell classification
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.bg`, color = snp.class)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd-16th HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd-16th HTO - colored by Souporcell classification (GEM1)") +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "snp.class")

# Compare UMI from each modality between Souporcell classes
Idents(gem1.hto) <- 'snp.class'
gem1.hto@active.ident <- factor(gem1.hto@active.ident, levels = c('singlet','multiplet','negative'))
VlnPlot(gem1.hto, features = 'nCount_HTO', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ADT', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_RNA', pt.size = 0.1, log = TRUE,  cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()

# 9) Analyze combined multiplet classification ----

# Create dataframe with combined metric multiplet classifications
pct.hto.df <- pct.hto.df %>%
  mutate(any.multip = case_when(
    hto.class == 'singlet' &
    snp.class == 'singlet' &
    scDblF.RNA.class == 'singlet' &
    scDblF.ATAC.class == 'singlet' ~ 'singlet',
    hto.class == 'multiplet' |
    snp.class == 'multiplet' |
    scDblF.RNA.class == 'multiplet' |
    scDblF.ATAC.class == 'multiplet' ~ 'multiplet',
    TRUE ~ 'negative'
  ))
pct.hto.df$any.multip <- factor(pct.hto.df$any.multip, levels = c('singlet','multiplet','negative'))

# Visualize results of combined multiplet assignment
pct.hto.df %>% ggplot(aes(x = any.multip, fill = any.multip)) + 
  geom_bar() + 
  scale_fill_manual(values = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = 'Combined multiplet assignment (GEM1)',
       x = 'Assignment',
       y = 'Number of cells') +
  theme_minimal()

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by combined multiplet assignment
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = any.multip)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO - colored by combined multiplet assignment (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "any.multip")

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by combined multiplet assignment
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = any.multip)) +
  geom_point(alpha = 0.25) +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO - colored by combined multiplet assignment (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "magenta", "singlet" = "lightgrey", "negative" = "green")) +
  coord_fixed() +
  labs(color = "any.multip")

# Visualize % UMI from maximum vs 2nd-16th HTO in each cell, colored by combined multiplet assignment
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.bg`, color = any.multip)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd-16th HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd-16th HTO - colored by combined multiplet assignment (GEM1)") +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "any.multip")

# Append combined multiplet assignments to Seurat object
pct.hto.df$any.multip <- factor(pct.hto.df$any.multip, levels = c('singlet','multiplet','negative'))
gem1.hto$any.multip <- pct.hto.df$any.multip

# Compare UMI from each modality between combined multiplet assignment classes
Idents(gem1.hto) <- 'any.multip'
gem1.hto@active.ident <- factor(gem1.hto@active.ident, levels = c('singlet','multiplet','negative'))
VlnPlot(gem1.hto, features = 'nCount_HTO', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ADT', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_RNA', pt.size = 0.1, log = TRUE,  cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()

# Find consensus multiplets between methods (expect few given differences between methods)
pct.hto.df <- pct.hto.df %>%
  mutate(joint.multip = case_when(
    hto.class == 'singlet' &
    snp.class == 'singlet' &
    scDblF.RNA.class == 'singlet' &
    scDblF.ATAC.class == 'singlet' ~ 'singlet',
    hto.class == 'multiplet' &
    snp.class == 'multiplet' &
    scDblF.RNA.class == 'multiplet' &
    scDblF.ATAC.class == 'multiplet' ~ 'multiplet',
    TRUE ~ 'negative'
  ))

# Append consensus multiplet assignments to Seurat object
pct.hto.df$joint.multip <- factor(pct.hto.df$joint.multip, levels = c('singlet','multiplet','negative'))
gem1.hto$joint.multip <- pct.hto.df$joint.multip

# Visualize % UMI from maximum vs secondary HTO in each cell, colored by consensus multiplet assignment
ggplot(pct.hto.df, aes(x = `pct.hto.max`, y = `pct.hto.2nd`, color = joint.multip)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(linewidth = 0.5, colour = 'black') +
  labs(x = "% UMI from Max HTO", y = "% UMI from 2nd HTO") +
  theme_minimal() +
  ggtitle("% UMI from Max vs 2nd HTO - colored by consensus multiplet assignment (GEM1)") +
  geom_hline(yintercept = 15) +
  geom_vline(xintercept = 50) +
  scale_color_manual(values = c("multiplet" = "red3", "singlet" = "skyblue2", "negative" = "grey")) +
  coord_fixed() +
  labs(color = "joint.multip")

# Compare UMI from each modality between consensus multiplet assignment classes
Idents(gem1.hto) <- 'joint.multip'
gem1.hto@active.ident <- factor(gem1.hto@active.ident, levels = c('singlet','multiplet','negative'))
VlnPlot(gem1.hto, features = 'nCount_HTO', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ADT', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_RNA', pt.size = 0.1, log = TRUE,  cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()
VlnPlot(gem1.hto, features = 'nCount_ATAC', pt.size = 0.1, log = TRUE, cols = c('singlet' = 'skyblue2', 'multiplet' = 'red3', 'negative' = 'grey')) + NoLegend()

# 10) Filter multiplets and unclassifiable cells from Seurat object ----

# Save Seurat object with all multiplet metrics appended before filtering
gem1.multip <- gem1.hto
saveRDS(gem1.multip, '/filepath/step8_filter_multiplets/gem1.multip.rds')

# Import multiplet object if starting from here
gem1.hto <- readRDS('/filepath/step8_filter_multiplets/gem1.multip.rds')

# Filter out multiplets and unclassified cells
gem1.singlet <- subset(gem1.hto, subset = 
                         hto.class == 'singlet' &
                         snp.class == 'singlet' &
                         scDblF.RNA.class == 'singlet' &
                         scDblF.ATAC.class == 'singlet'
)

# 11) Filter possibly missed multiplets with exceptionally high UMI counts ----

# Survey 99th percentile of UMI counts across modalities within filtered 'singlets'
quantile(gem1.singlet$nCount_HTO, 0.999)
quantile(gem1.singlet$nCount_ADT, 0.999)
quantile(gem1.singlet$nCount_RNA, 0.999)
quantile(gem1.singlet$nFeature_RNA, 0.999)
quantile(gem1.singlet$nCount_ATAC, 0.999)
quantile(gem1.singlet$nFeature_ATAC, 0.999)

# Visualize outlier cutoff for UMI counts across modalities
Idents(gem1.singlet) <- 'hto.sort'
VlnPlot(gem1.singlet, features = 'nCount_HTO', pt.size = 0.1, log = TRUE) + NoLegend() + geom_hline(yintercept = 4000)
VlnPlot(gem1.singlet, features = 'nCount_ADT', pt.size = 0.1, log = TRUE) + NoLegend() + geom_hline(yintercept = 4500)
VlnPlot(gem1.singlet, features = 'nCount_RNA', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 20000)
VlnPlot(gem1.singlet, features = 'nFeature_RNA', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 5000)
VlnPlot(gem1.singlet, features = 'nCount_ATAC', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 50000)
VlnPlot(gem1.singlet, features = 'nFeature_ATAC', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 15000)

# Filter outlier UMI droplets from object
gem1.singlet.umi <- subset(gem1.singlet, subset = 
                             nCount_HTO < 4000 &
                             nCount_ADT < 4500 &
                             nCount_RNA < 20000 &
                             nFeature_RNA < 5000 &
                             nCount_ATAC < 50000 &
                             nFeature_ATAC < 15000)

# Determine number of multiplets and outlier UMI droplets removed
(ncol(gem1.hto)) # 5120 input
(ncol(gem1.singlet)) # 4033 singlets
(ncol(gem1.singlet) - ncol(gem1.singlet.umi)) # 13 removed
(ncol(gem1.singlet.umi)) # 4020 remaining

# Visualize new UMI distribution across modalities
VlnPlot(gem1.singlet.umi, features = 'nCount_HTO', pt.size = 0.1, log = TRUE) + NoLegend() + geom_hline(yintercept = 4000) + ylim(0,4500)
VlnPlot(gem1.singlet.umi, features = 'nCount_ADT', pt.size = 0.1, log = TRUE) + NoLegend() + geom_hline(yintercept = 4500) + ylim(0,5000)
VlnPlot(gem1.singlet.umi, features = 'nCount_RNA', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 20000) + ylim(0,20500)
VlnPlot(gem1.singlet.umi, features = 'nFeature_RNA', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 5000) + ylim(0,5500)
VlnPlot(gem1.singlet.umi, features = 'nCount_ATAC', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 50000) + ylim(0,50500)
VlnPlot(gem1.singlet.umi, features = 'nFeature_ATAC', pt.size = 0.1, log = FALSE) + NoLegend() + geom_hline(yintercept = 15000) + ylim(0,15500)

# 12) Save singlet Seurat object ----
saveRDS(gem1.singlet.umi, '/filepath/step8_filter_multiplets/gem1.singlet.umi.rds')

# 13) Repeat for all 8 GEM wells and proceed to next step ----