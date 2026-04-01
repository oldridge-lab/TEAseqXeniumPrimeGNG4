# Xenium Prime Data Preprocessing Step 1
# For each sample, import data from Xenium Onboard Analysis and create Seurat object

# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 
# Sample - FFPE Tonsil, Donor TC656A, 4YO M with OSA (Table S1, six total samples)

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(arrow) # 18.1.0.1

# Set working directory
setwd("/filepath/xenium_data_processing/step1_preprocess/")

# Define path to outputs from Xenium Onboard Analysis for selected sample
tc656a_path <- "/filepath/xenium_data_processing/xenium_xoa_output_folder/tc656a/"

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
tc656a.data <- ReadXenium(tc656a_path, outs = c('matrix','microns'),
                          type = c('segmentations', 'centroids'),
                          mols.qv.threshold = 20)

# Finish Seurat object creation using steps normally in LoadXenium function
tc656a.segmentations <- list(centroids = CreateCentroids(tc656a.data$centroids), 
                             segmentation = CreateSegmentation(tc656a.data$segmentations))
tc656a.coords <- CreateFOV(coords = tc656a.segmentations, 
                           type = c("segmentation","centroids"),
                           molecules = tc656a.data$microns, assay = "RNA")
tc656a.obj <- CreateSeuratObject(counts = tc656a.data$matrix[["Gene Expression"]], assay = "RNA")
tc656a.obj[["fov"]] <- tc656a.coords

# Add quality control data from XOA as separate assays to object
tc656a.obj[["BlankCodeword"]] <- CreateAssayObject(counts = tc656a.data$matrix[["Unassigned Codeword"]])
tc656a.obj[["ControlCodeword"]] <- CreateAssayObject(counts = tc656a.data$matrix[["Negative Control Codeword"]])
tc656a.obj[["ControlProbe"]] <- CreateAssayObject(counts = tc656a.data$matrix[["Negative Control Probe"]])
tc656a.obj[["GenomicControl"]] <- CreateAssayObject(counts = tc656a.data$matrix[["Genomic Control"]])
tc656a.obj[["DeprecatedCodeword"]] <- CreateAssayObject(counts = tc656a.data$matrix[["Deprecated Codeword"]])

# Add donor metadata
tc656a.obj@project.name <- 'TC656A'
tc656a.obj$donor <- 'TC656A'
tc656a.obj$age <- 4
tc656a.obj$asab <- 'Male'
tc656a.obj$indication <- 'OSA'
tc656a.obj$slide <- 'Slide1'
Idents(tc656a.obj) <- 'donor'

# Save preprocessed object 
saveRDS(tc656a.obj, file = '/filepath/xenium_data_processing/step1_preprocess/tc656a.preproc.seurat.obj.rds')

# Repeat for all six samples