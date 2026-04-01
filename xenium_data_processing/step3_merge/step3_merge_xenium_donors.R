# Xenium Prime Data Preprocessing Step 3
# Merging Seurat objects from each Xenium Prime tonsil sample (Table S1) after quality control filtering
# Assay - 5001-plex Xenium Prime Spatial Transcriptomics with 100-plex custom add-on panel (Table S9) and Multimodal Segmentation 

# Set up R working environment

# Set seed
set.seed(26)

# Load packages and record version numbers
library(Seurat) # 5.1.0

# Set working directory
setwd("/filepath/xenium_data_processing/step3_merge/")

# Import Seurat objects for each of the six tonsil FFPE Xenium Prime samples after quality control filtering in Step 2
print('importing objects')
tc656a.qc <- readRDS(file = '/filepath/xenium_data_processing/step2_qc/tc656a_post_qc_obj.rds')
tc654b.qc <- readRDS(file = '/filepath/xenium_data_processing/step2_qc/tc654b_post_qc_obj.rds')
tc653b.qc <- readRDS(file = '/filepath/xenium_data_processing/step2_qc/tc653b_post_qc_obj.rds')
tc661b.qc <- readRDS(file = '/filepath/xenium_data_processing/step2_qc/tc661b_post_qc_obj.rds')
tc659a.qc <- readRDS(file = '/filepath/xenium_data_processing/step2_qc/tc659a_post_qc_obj.rds')
tc657b.qc <- readRDS(file = '/filepath/xenium_data_processing/step2_qc/tc657b_post_qc_obj.rds')
print('all objects imported')

# Verify cumulative cell number
cat((ncol(tc656a.qc) + ncol(tc654b.qc) + ncol(tc653b.qc) + ncol(tc661b.qc) + ncol(tc659a.qc) + ncol(tc657b.qc)),'cells total expected in post-QC merged object')

# Merge objects, appending sample identifier to cell barcodes
print('merging objects')
xp.merged.obj <- merge(tc656a.qc, 
                        y = c(tc654b.qc, 
                              tc653b.qc,
                              tc661b.qc,
                              tc659a.qc,
                              tc657b.qc),
                        add.cell.ids = c("TC656A",
                                         "TC654B",
                                         "TC653B",
                                         "TC661B",
                                         "TC659A",
                                         "TC657B"),
                        project = "XP_Merged")
print('objects merged')

# Save merged object
print('saving post-QC merged object')
saveRDS(xp.merged.obj, '/filepath/xenium_data_processing/step3_merge/xp_obj_step3.rds')
print('merged object saved')

print('R script done')