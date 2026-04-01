# TEAseq Preprocessing Step 5 (GEM7) ----

# Substeps 
# 1) Set up R working environment
# 2) Import cellranger-arc filtered feature barcode matrix and cellranger multi raw feature barcode matrix
# 3) Remove ambient RNA using SoupX
# 4) Create Seurat object to which we will add all other modalities (ATAC, ADT, HTO, and SNP)
# 5) Create Signac ATAC assay and to Seurat object
# 6) Subset object on filtered RNA/ATAC cell barcodes that are found in raw ADT/HTO feature matrix (should contain 100% overlap)
# 7) Add HTO data to Seurat object
# 8) Add ADT data to Seurat object
# 9) Add SNP-based demultiplexing metadata from Souporcell to Seurat object
# 10) Save pre-processed Seurat object
# 11) Repeat for all GEM wells, then proceed to next step

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

# Set working directory to save objects in pre-processing folder
setwd("/filepath/step5_preprocess_gems")

# 2) Import cell x feature count matrices from cellranger-arc and cellranger multi ----

# Import cellranger multi raw RNA and ADT/HTO feature barcode matrix
multi_7.data <- Read10X("/filepath/cellranger_multi_output_newlib/GEM7_multi/outs/multi/count/raw_feature_bc_matrix", unique.features = TRUE)

# Import cellranger-arc count filtered RNA and ATAC feature barcode matrix
arc_7.data <- Read10X("/filepath/cellranger-arc_output/GEM7_ARC/outs/filtered_feature_bc_matrix", unique.features = TRUE)

# 3) Remove ambient RNA using SoupX----

# Load in cellranger-arc count filtered count matrix in SoupX-readable format
soup_7 <- load10X("/filepath/cellranger-arc_output/GEM7_ARC/outs/", use.names = TRUE, unique.features = TRUE)

# Perform automatic estimation of ambient RNA
soup_7 <- autoEstCont(soup_7)

# Adjust RNA count matrix to remove estimated ambient counts
clean_7 <- adjustCounts(soup_7)

# 4) Create Seurat object for GEM7 using cleaned RNA matrix, to which we will add all other modalities ----
obj_gem7 <- CreateSeuratObject(counts = clean_7, project = "GEM_7")

# 5) Create Signac ATAC assay and add to Seurat object ----

# Access cellranger-arc filtered cell barcode x Tn5 integration peak matrix
atac_7.counts <- arc_7.data$Peaks

# Record path to fragment file (full list of unique ATAC fragments per single cell barcode)
frag_7.file <- "/filepath/cellranger-arc_output/GEM7_ARC/outs/atac_fragments.tsv.gz"

# Retrieve gene annotations for Genome Reference Consortium Human Build 38 (GRCh38)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# Create Signac ATAC assay using cell x peak matrix, fragment list, and GRCh38 annotations
chrom_assay_7 <- CreateChromatinAssay(
  counts = atac_7.counts,
  sep = c(":", "-"),
  fragments = frag_7.file,
  annotation = annotations
)

# Add ATAC assay to Seurat object
obj_gem7[["ATAC"]] <- chrom_assay_7

# 6) Subset on cell barcodes shared between filtered ATAC/RNA feature matrix (cellranger-arc output) and raw ADT/HTO feature matrix (cellranger multi output) ----

# Extract barcodes from each cellranger output feature matrix
multi_7_barcodes <- colnames(multi_7.data[["Gene Expression"]])
arc_7_barcodes <- colnames(arc_7.data[["Gene Expression"]])

# Find barcode intersect between cellranger-arc and cellranger multi output matrices
common_7_barcodes <- intersect(multi_7_barcodes, arc_7_barcodes)

# Subset Seurat object on common barcodes
obj_gem7_subset <- subset(obj_gem7, cells = common_7_barcodes)
obj_gem7_subset$GEM <- "GEM_7"

# 7) Add HTO assay to Seurat object ----

# Create cell x HTO feature matrix using cellranger multi output
hto_7_count <- multi_7.data
hto_7_count[["Antibody Capture"]] <- hto_7_count[["Antibody Capture"]][grep("HTO", rownames(hto_7_count[["Antibody Capture"]])), ]
hto_7_count[["Gene Expression"]] <- NULL
hto_7 <- CreateSeuratObject(counts = hto_7_count)
hto_7_subset <- subset(hto_7, cells = common_7_barcodes)
RenameAssays(hto_7_subset, assay.name = 'RNA', new.assay.name = 'HTO')
obj_gem7_subset[["HTO"]] <- CreateAssayObject(counts = hto_7_subset@assays$RNA$counts)  

# 8) Add ADT assay to Seurat object ----

# Create cell x ADT feature matrix using cellranger multi output
adt_7_count <- multi_7.data
adt_7_count[["Antibody Capture"]] <- adt_7_count[["Antibody Capture"]][-grep("HTO", rownames(adt_7_count[["Antibody Capture"]])), ]
adt_7_count[["Gene Expression"]] <- NULL
adt_7 <- CreateSeuratObject(counts = adt_7_count)
adt_7_subset <- subset(adt_7, cells = common_7_barcodes)
RenameAssays(adt_7_subset, assay.name = 'RNA', new.assay.name = 'ADT')
obj_gem7_subset[["ADT"]] <- CreateAssayObject(counts = adt_7_subset@assays$RNA$counts)  

# 9) Add SNP-based demultiplexing metadata from Souporcell to Seurat object ----

# Import souporcell output SNP-based sample demultiplexing cluster file
snp.metadata.gem7 <- as.data.frame(read_tsv("/filepath/step4_snp_demultiplex/gem7/clusters.tsv"))

# Add souporcell output metadata to Seurat object
obj_gem7_subset <- AddMetaData(obj_gem7_subset, metadata = snp.metadata.gem7)

# 10) Save pre-processed GEM object ----
saveRDS(obj_gem7_subset, "/filepath/step5_preprocess_gems/gem7_preprocessed.rds")

# 11) Repeat for all GEM wells, then proceed to Step 6 ----