# Xenium Prime Data Preprocessing Step 10
# Annotation of cell neighborhoods generated in Step 9
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

# Cell neighborhood (CN) analysis used L1 cluster annotations as input (not L2 T cell subclusters), k = 40 neighbors, and n = 10 neighborhoods
# Original method from Schürch, Christian M. et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. Cell, Volume 182, Issue 5, 1341 - 1359.e19
# Code adapted from related GitHub page - https://github.com/nolanlab/NeighborhoodCoordination
# CN generation for tonsil Xenium Prime dataset was performed in collaboration with Yutong Zhu from Oldridge Lab. Refer to Python scripts in 'step9_wsi_registration_cn_analysis' folder for details regarding CN generation.

# CNs were annotated as below, then visualized using Xenium Explorer. Refer to Fig S16 and related script for further data supporting annotation of CNs.

# Set up R working environment ----

# Set seed
set.seed(26)

# Set working directory
setwd("/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/")

# Load packages and record version numbers
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(BPCells) # 0.3.0
library(presto) # 1.0.0
library(tidyverse) # 2.0.0

# Import L1-3 Seurat objects from Step 8C annotation ----

# L1 object
xp.obj <- readRDS(file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/xp_l1_obj_step8c.rds')

# L2 nnCD4 T cell subclustering object
xp_nncd4t_obj <- readRDS(file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/xp_l2_nncd4t_obj_step8c.rds')

# L3 Tfh-only subset object
xp_tfh_obj <- readRDS(file = '/filepath/xenium_data_processing/step8_l2_nncd4t_subclustering/xp_l3_tfh_subset_obj_step8c.rds')

# Import original CN and closed GC CN data from Step 9 ----

# Import original CN assignment for each segmented cell barcode including sample, L1 cluster, and cluster_id metadata
cell_neighb_csv <- read.csv('/filepath/cell_neighborhood_analysis_output_folder/neighborhoods_k_40.csv')

# Secondary 'GC Closed' CN analysis, wherein CN2/5 LZ/DZ were merged and any 'holes' enclosed by the merged GC CNs are filled 
cn_gc_closed_csv <- read.csv('/filepath/closed_gc_cn_analysis_output_folder/neighborhoods_k_40_closed.csv')

# Annotate original numerical CNs - CN0-CN9 ----

# Save numerical CN assignments as metadata column in L1 XP Seurat object
xp.obj$cn <- cell_neighb_csv$neighborhoodk40n10

# Annotate numerical CNs
cn_annot <- c(
  '0'  = 'CN0 Fol Border',
  '1'  = 'CN1 TCZ Inner',
  '2'  = 'CN2 GC LZ',
  '3'  = 'CN3 TCZ Outer',
  '4'  = 'CN4 Epi Outer ',
  '5'  = 'CN5 GC DZ',
  '6'  = 'CN6 Mantle',
  '7'  = 'CN7 Epi Inner', 
  '8'  = 'CN8 ASC Rich',
  '9'  = 'CN9 Connective'
)

# Transfer CN annotations to L1 object
xp.obj$cn_annot <- xp.obj$cn
Idents(xp.obj) <- 'cn_annot'
xp.obj <- RenameIdents(xp.obj, cn_annot)
xp.obj$cn_annot <- Idents(xp.obj)

# Transfer CN annotations to L2 nnCD4 T cell object
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_nncd4t_obj$cn <- xp.obj@meta.data[colnames(xp_nncd4t_obj), "cn", drop = TRUE]
xp_nncd4t_obj$cn_annot <- xp.obj@meta.data[colnames(xp_nncd4t_obj), "cn_annot", drop = TRUE]

# Transfer CN annotations to L3 Tfh-only subset object
DefaultAssay(xp_tfh_obj) <- 'RNA'
xp_tfh_obj$cn <- xp.obj@meta.data[colnames(xp_tfh_obj), "cn", drop = TRUE]
xp_tfh_obj$cn_annot <- xp.obj@meta.data[colnames(xp_tfh_obj), "cn_annot", drop = TRUE]

# Annotate CNs from secondary merged GC CN analysis ----

# Save numerical CN assignments as metadata column in L1 XP Seurat object
xp.obj$cn_gc_closed <- cn_gc_closed_csv$neighborhoodk40n10_gcClosed
table(xp.obj$cn_gc_closed)

# Annotation of CN_GC_Closed
cn_gc_closed_annot <- c(
  '0'  = 'CN0 Fol Border',
  '1'  = 'CN1 TCZ Inner',
  '2'  = 'CN2 GC', # Compared to original CN annotation, CN2 LZ and CN5 DZ now merged into GC CN, along with any GC-enclosed cells initially assigned to other CN 
  '3'  = 'CN3 TCZ Outer',
  '4'  = 'CN4 Epi Outer ',
  '6'  = 'CN6 Mantle',
  '7'  = 'CN7 Epi Inner', 
  '8'  = 'CN8 ASC Rich',
  '9'  = 'CN9 Connective'
)

# Transfer closed GC CN annotations to L1 object
xp.obj$cn_gc_closed_annot <- xp.obj$cn_gc_closed
Idents(xp.obj) <- 'cn_gc_closed_annot'
xp.obj <- RenameIdents(xp.obj, cn_gc_closed_annot)
xp.obj$cn_gc_closed_annot <- Idents(xp.obj)
table(xp.obj$cn_gc_closed_annot)

# Transfer closed GC CN annotations to L2 nnCD4 T cell subclustering object
DefaultAssay(xp_nncd4t_obj) <- 'RNA'
xp_nncd4t_obj$cn_gc_closed <- xp.obj@meta.data[colnames(xp_nncd4t_obj), "cn_gc_closed", drop = TRUE]
xp_nncd4t_obj$cn_gc_closed_annot <- xp.obj@meta.data[colnames(xp_nncd4t_obj), "cn_gc_closed_annot", drop = TRUE]
table(xp_nncd4t_obj$cn_gc_closed_annot)

# Transfer closed GC CN annotations to L3 Tfh-only subset object
DefaultAssay(xp_tfh_obj) <- 'RNA'
xp_tfh_obj$cn_gc_closed <- xp.obj@meta.data[colnames(xp_tfh_obj), "cn_gc_closed", drop = TRUE]
xp_tfh_obj$cn_gc_closed_annot <- xp.obj@meta.data[colnames(xp_tfh_obj), "cn_gc_closed_annot", drop = TRUE]
table(xp_tfh_obj$cn_gc_closed_annot)

# Export numerical CN assignments per cell for each sample and export to explore using Xenium Explorer ----

# TC653B - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_number_per_sample_csv/tc653b_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC656A - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn'
cells_tc656a <- colnames(xp.obj)[ xp.obj$donor == "TC656A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc656a)
idents_tc656a <- Idents(xp.obj)[ cells_tc656a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc656a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_number_per_sample_csv/tc656a_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC654B - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn'
cells_tc654b <- colnames(xp.obj)[ xp.obj$donor == "TC654B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc654b)
idents_tc654b <- Idents(xp.obj)[ cells_tc654b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc654b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_number_per_sample_csv/tc654b_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC661B - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn'
cells_tc661b <- colnames(xp.obj)[ xp.obj$donor == "TC661B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc661b)
idents_tc661b <- Idents(xp.obj)[ cells_tc661b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc661b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_number_per_sample_csv/tc661b_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC659A - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn'
cells_tc659a <- colnames(xp.obj)[ xp.obj$donor == "TC659A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc659a)
idents_tc659a <- Idents(xp.obj)[ cells_tc659a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc659a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_number_per_sample_csv/tc659a_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC657B - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn'
cells_tc657b <- colnames(xp.obj)[ xp.obj$donor == "TC657B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc657b)
idents_tc657b <- Idents(xp.obj)[ cells_tc657b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc657b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_number_per_sample_csv/tc657b_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# Export annotated CN assignments per cell for each sample and export to explore using Xenium Explorer ----

# TC653B - extract CN annotation and save as csv to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_annot'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_annot_per_sample_csv/tc653b_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC656A - extract CN annotation and save as csv to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_annot'
cells_tc656a <- colnames(xp.obj)[ xp.obj$donor == "TC656A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc656a)
idents_tc656a <- Idents(xp.obj)[ cells_tc656a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc656a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_annot_per_sample_csv/tc656a_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC654B - extract CN annotation and save as csv to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_annot'
cells_tc654b <- colnames(xp.obj)[ xp.obj$donor == "TC654B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc654b)
idents_tc654b <- Idents(xp.obj)[ cells_tc654b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc654b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_annot_per_sample_csv/tc654b_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC661B - extract CN annotation and save as csv to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_annot'
cells_tc661b <- colnames(xp.obj)[ xp.obj$donor == "TC661B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc661b)
idents_tc661b <- Idents(xp.obj)[ cells_tc661b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc661b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_annot_per_sample_csv/tc661b_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC659A - extract CN annotation and save as csv to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_annot'
cells_tc659a <- colnames(xp.obj)[ xp.obj$donor == "TC659A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc659a)
idents_tc659a <- Idents(xp.obj)[ cells_tc659a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc659a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_annot_per_sample_csv/tc659a_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC657B - extract CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_annot'
cells_tc657b <- colnames(xp.obj)[ xp.obj$donor == "TC657B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc657b)
idents_tc657b <- Idents(xp.obj)[ cells_tc657b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc657b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/cn_annot_per_sample_csv/tc657b_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)


# Export numerical closed GC CN assignments per cell for each sample to explore in Xenium Explorer ----

# TC653B - extract closed GC CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_number_per_sample_csv/tc653b_gc_closed_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC656A - extract closed GC CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed'
cells_tc656a <- colnames(xp.obj)[ xp.obj$donor == "TC656A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc656a)
idents_tc656a <- Idents(xp.obj)[ cells_tc656a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc656a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_number_per_sample_csv/tc656a_gc_closed_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC654B - extract closed GC CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed'
cells_tc654b <- colnames(xp.obj)[ xp.obj$donor == "TC654B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc654b)
idents_tc654b <- Idents(xp.obj)[ cells_tc654b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc654b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_number_per_sample_csv/tc654b_gc_closed_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC661B - extract closed GC CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed'
cells_tc661b <- colnames(xp.obj)[ xp.obj$donor == "TC661B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc661b)
idents_tc661b <- Idents(xp.obj)[ cells_tc661b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc661b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_number_per_sample_csv/tc661b_gc_closed_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC659A - extract closed GC CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed'
cells_tc659a <- colnames(xp.obj)[ xp.obj$donor == "TC659A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc659a)
idents_tc659a <- Idents(xp.obj)[ cells_tc659a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc659a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_number_per_sample_csv/tc659a_gc_closed_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC657B - extract closed GC CN number and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed'
cells_tc657b <- colnames(xp.obj)[ xp.obj$donor == "TC657B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc657b)
idents_tc657b <- Idents(xp.obj)[ cells_tc657b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc657b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_number_per_sample_csv/tc657b_gc_closed_cn_numb_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# Export closed GC CN annotations per cell for each sample to explore in Xenium Explorer ----

# TC653B - extract closed GC CN annotation and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed_annot'
cells_tc653b <- colnames(xp.obj)[ xp.obj$donor == "TC653B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc653b)
idents_tc653b <- Idents(xp.obj)[ cells_tc653b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc653b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_annot_per_sample_csv/tc653b_gc_closed_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC656A - extract closed GC CN annotation and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed_annot'
cells_tc656a <- colnames(xp.obj)[ xp.obj$donor == "TC656A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc656a)
idents_tc656a <- Idents(xp.obj)[ cells_tc656a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc656a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_annot_per_sample_csv/tc656a_gc_closed_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC654B - extract closed GC CN annotation and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed_annot'
cells_tc654b <- colnames(xp.obj)[ xp.obj$donor == "TC654B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc654b)
idents_tc654b <- Idents(xp.obj)[ cells_tc654b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc654b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_annot_per_sample_csv/tc654b_gc_closed_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC661B - extract closed GC CN annotation and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed_annot'
cells_tc661b <- colnames(xp.obj)[ xp.obj$donor == "TC661B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc661b)
idents_tc661b <- Idents(xp.obj)[ cells_tc661b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc661b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_annot_per_sample_csv/tc661b_gc_closed_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC659A - extract closed GC CN annotation and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed_annot'
cells_tc659a <- colnames(xp.obj)[ xp.obj$donor == "TC659A" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc659a)
idents_tc659a <- Idents(xp.obj)[ cells_tc659a ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc659a),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_annot_per_sample_csv/tc659a_gc_closed_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# TC657B - extract closed GC CN annotation and save as CSV to import into Xenium Explorer
DefaultAssay(xp.obj) <- 'RNA'
Idents(xp.obj) <- 'cn_gc_closed_annot'
cells_tc657b <- colnames(xp.obj)[ xp.obj$donor == "TC657B" ]
trimmed_ids <- sub("^[^_]+_", "", cells_tc657b)
idents_tc657b <- Idents(xp.obj)[ cells_tc657b ]
write.csv(
  data.frame(
    cell_id = trimmed_ids,
    group   = as.character(idents_tc657b),
    stringsAsFactors = FALSE
  ),
  file      = "/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/gc_closed_cn_annot_per_sample_csv/tc657b_gc_closed_cn_annot_n10_k40.csv",
  row.names = FALSE,
  quote     = FALSE
)

# Save Seurat objects with CN information for both original and 'GC closed' analyses  ----

# L1 object
saveRDS(xp.obj, file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l1_obj_step10.rds')

# L2 nnCD4 T subclustering object 
saveRDS(xp_nncd4t_obj, file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l2_nncd4t_obj_step10.rds')

# L3 Tfh-only subset object
saveRDS(xp_tfh_obj, file = '/filepath/xenium_data_processing/step10_cell_neighborhood_annotation/xp_l3_tfh_subset_obj_step10.rds')