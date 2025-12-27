# Xenium Prime Data Preprocessing Step 1
# For each sample, import data from Xenium Onboard Analysis and create Seurat object

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC659A, 4YO M with SDB (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(arrow) # 18.1.0.1

# Set working directory
setwd("/filepath/xenium_data_processing/step1_preprocess/")

# Define path to outputs from Xenium Onboard Analysis for selected sample
tc659a_path <- "/filepath/xenium_data_processing/xenium_xoa_output_folder/tc659a/"

# Define function to import data from XOA
# LoadXenium function was not compatible with XOA v3.0 upon release - fixed in Seurat v5.2.0 (e.g. obj <- LoadXenium(path, fov = "fov"))
# For Xenium data analysis we maintained using Seurat v5.1.0 consistent with our TEAseq analysis environment for this project (https://github.com/satijalab/seurat/releases)
# Therefore, we defined a XOA data reading function that could accommodate Seurat v5.1.0 based on community feedback  (https://github.com/satijalab/seurat/issues/9060)

ReadXenium <- function (data.dir, outs = c("matrix", "microns"), type = "centroids", 
                        mols.qv.threshold = 20) 
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location, 
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}

# Use new read function to import Xenium data
tc659a.data <- ReadXenium(tc659a_path, outs = c('matrix','microns'),
                          type = c('segmentations', 'centroids'),
                          mols.qv.threshold = 20)

# Finish Seurat object creation using steps normally in LoadXenium function
tc659a.segmentations <- list(centroids = CreateCentroids(tc659a.data$centroids), 
                             segmentation = CreateSegmentation(tc659a.data$segmentations))
tc659a.coords <- CreateFOV(coords = tc659a.segmentations, 
                           type = c("segmentation","centroids"),
                           molecules = tc659a.data$microns, assay = "RNA")
tc659a.obj <- CreateSeuratObject(counts = tc659a.data$matrix[["Gene Expression"]], assay = "RNA")
tc659a.obj[["fov"]] <- tc659a.coords

# Add quality control data from XOA as separate assays to object
tc659a.obj[["BlankCodeword"]] <- CreateAssayObject(counts = tc659a.data$matrix[["Unassigned Codeword"]])
tc659a.obj[["ControlCodeword"]] <- CreateAssayObject(counts = tc659a.data$matrix[["Negative Control Codeword"]])
tc659a.obj[["ControlProbe"]] <- CreateAssayObject(counts = tc659a.data$matrix[["Negative Control Probe"]])
tc659a.obj[["GenomicControl"]] <- CreateAssayObject(counts = tc659a.data$matrix[["Genomic Control"]])
tc659a.obj[["DeprecatedCodeword"]] <- CreateAssayObject(counts = tc659a.data$matrix[["Deprecated Codeword"]])

# Add donor metadata
tc659a.obj@project.name <- 'TC659A'
tc659a.obj$donor <- 'TC659A'
tc659a.obj$age <- 4
tc659a.obj$asab <- 'Male'
tc659a.obj$indication <- 'SDB'
tc659a.obj$slide <- 'Slide2'
Idents(tc659a.obj) <- 'donor'

# Save preprocessed object 
saveRDS(tc659a.obj, file = '/filepath/xenium_data_processing/step1_preprocess/tc659a.preproc.seurat.obj.rds')

# Repeat for all six samples