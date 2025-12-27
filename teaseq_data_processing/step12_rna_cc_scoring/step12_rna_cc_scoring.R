# Step 12 - Cell-cycle scoring of merged object (RNA modality) ----

# Substeps
# 1) Set up R working environment
# 2) Import 8 x GEM merged Seurat object
# 3) Perform cell-cycle gene scoring
# 4) Calculate difference between S and G2M scores
# 5) Rename cell-cycle phase metadata column
# 6) Record cell-cycle scoring results
# 7) Save merged Seurat object with cell-cycle scores and phase metadata appended
# 8) Proceed to next step - Level 1 TEAseq object dimensionality reduction and 3WNN clustering ----

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

# Set working directory to save objects created during this substep
setwd("/filepath/step12_rna_cc_scoring")

# 2) Import 8 x GEM merged Seurat object  ----

merged.obj <- readRDS("/filepath/step11_merge_gems/merged.obj.rds")

# 3) Perform cell-cycle gene scoring ----

DefaultAssay(merged.obj) <- "RNA"
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
merged.obj <- NormalizeData(merged.obj)
merged.cc.scored <- CellCycleScoring(merged.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 4) Calculate difference between S and G2M Phase scores ----

# Can regress CC.Difference in future steps to preserve differences between cycling and non-cycling cells, while regressing differences among cycling cells. Ultimately no regression of cycling signature was performed.
# Refer to https://satijalab.org/seurat/articles/cell_cycle_vignette#alternate-workflow

merged.cc.scored$CC.Difference <- merged.cc.scored$S.Score - merged.cc.scored$G2M.Score

# 5) Rename cell-cycle phase metadata column to be more specific ----
merged.cc.scored$CC.Phase <- merged.cc.scored$Phase
merged.cc.scored$Phase <- NULL # removing old default metadata column name

# 6) Record cell-cycle scoring results ----
table(merged.cc.scored$CC.Phase)

# G1 - 11565 
# G2M - 11883  
# S - 8758

# 7) Save merged Seurat object with cell-cycle scores and phase metadata appended ----
saveRDS(merged.cc.scored, "/filepath/step12_rna_cc_scoring/merged.cc.scored.rds")

# 8) Proceed to next step - L1 clustering and 3WNN analysis ----