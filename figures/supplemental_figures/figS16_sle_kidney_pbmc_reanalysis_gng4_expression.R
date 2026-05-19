# Sam Barnett Dubensky et al.
# Derek A. Oldridge & Laura A. Vella Labs at the Children's Hospital of Philadelphia
# Multimodal analysis defines GNG4 as a distinguishing feature of germinal center-positioned Tfh in human lymphoid tissue
# Code and data visualization for Fig. S16C-S16E (related to Fig. 6D-6F)
# Fig. S16 - GNG4 risk variant mapping in rheumatoid arthritis and expression by cTfh-like cells in systemic lupus erythematosus, related to Figure 6

# Set up R working environment ----

# Set working directory to Fig 6
setwd('/filepath/fig6/fig6_supp/')

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
library(ggseqlogo) # 0.2
library(caret) # 6.0-94
library(presto) # 1.0.0
library(SoupX) # 1.6.2
library(scDblFinder) # 1.19.1
library(Nebulosa) # 1.14.0
library(clusterProfiler) # 4.12.6 
library(enrichplot) # 1.24.2
library(org.Hs.eg.db) # 3.19.1
library(harmony) # 1.2.0
library(SCpubr) # 2.0.2
library(paletteer) # 1.6.0
library(RColorBrewer) # 1.1-3
library(EnhancedVolcano) # 1.22.0
library(AUCell) # 1.26.0
library(msigdbr) # 7.5.1
library(RcisTarget) # 1.23.1
library(GENIE3) # 1.26.0
library(SCENIC) # 1.3.1
library(arrow) # 18.1.0.1
library(doMC) # 1.3.8
library(R2HTML) # 2.3.4
library(reshape2) # 1.4.4
library(ggtern) # 3.5.0
library(ggh4x) # 0.3.0
library(SeuratWrappers) # 0.3.5
library(monocle3) # 1.3.7
library(cisTopic) # 0.3.0
library(densityClust) # 0.3.3
library(Rtsne) # 0.17
library(scatterplot3d) # 1.2
library(fastcluster) # 1.2.6
library(rtracklayer) # 1.64.0
library(tidyr) # 1.3.1
library(tidyverse) # 2.0.0
library(pheatmap) # 1.0.12
library(ggraph) # 2.2.1
library(igraph) # 2.0.3
library(BPCells) # 0.3.0
library(future) # 1.33.2
library(FNN) # 1.1.4
library(tibble) # 3.2.1
library(purrr) # 1.0.2
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # 3.18.0
library(ComplexHeatmap) # 2.20.0
library(circlize) # 0.4.16
library(BiocParallel) # 1.38.0
library(writexl) # 1.5.4
library(ggnewscale) # 0.5.0
library(readxl) # 1.4.3
library(UCell) # 2.8.0
library(scales) # 1.3.0
library(rlang) # 1.1.4
library(ggrepel) # 0.9.5
library(speckle) # 1.4.0

# Fig. S16C - recreate Seurat object from Al Souz et al.

# Reanalysis of kidney tissue scRNAseq data from adult healthy donors versus systemic lupus erythematosus (SLE) patients 
# Originally published by Al Souz et al. (2026) Immunity
# https://pubmed.ncbi.nlm.nih.gov/42013863/
# Organ injury in systemic autoimmunity is mediated by stem-like CD8+ T cells arising from tissue-draining lymph nodes

# Retrieve scRNAseq files and metadata from Broad Single Cell Portal - https://singlecell.broadinstitute.org/single_cell/study/SCP3488/amp2-lupus-kidney-single-cell

# Define path to folders containing Al Souz et al. files
base_dir <- "/filepath/fig6/SCP3488/expression/69c1753e771a5b0296040b35/single_cell_portal/expr_raw"

# Assemble cell-by-gene count matrix
counts <- ReadMtx(
  mtx = file.path(base_dir, "matrix_raw.mtx.gz"),
  features = file.path(base_dir, "features_raw.tsv.gz"),
  cells = file.path(base_dir, "barcodes_raw.tsv.gz"),
  feature.column = 1
)

# Create Seurat object
kidney_sle_obj <- CreateSeuratObject(
  counts = counts,
  project = "SCP3488"
)

# Import metadata from Al Souz et al.
md <- fread("/filepath/fig6/SCP3488/metadata/single_cell_portal/metadata.txt.gz", data.table = FALSE)

# Match cell barcode rows between metadata and Seurat object
rownames(md) <- md$NAME
md <- md[colnames(kidney_sle_obj), , drop = FALSE]

# Add metadata to Seurat object
kidney_sle_obj <- AddMetaData(kidney_sle_obj, metadata = md)

# Import UMAP reduction coordinates from Al Souz et al.
umap_tnk <- fread("/filepath/fig6/SCP3488/cluster/single_cell_portal/umap_TNK.txt.gz", data.table = FALSE)
colnames(umap_tnk)
head(umap_tnk)

# Match UMAP dataframe cell barcodes to Seurat object
rownames(umap_tnk) <- umap_tnk$NAME
umap_df <- umap_tnk[colnames(kidney_sle_obj), c('X', 'Y'), drop = FALSE]
umap_mat <- as.matrix(sapply(umap_df, as.numeric))
rownames(umap_mat) <- colnames(kidney_sle_obj)
colnames(umap_mat) <- c("UMAPTNK_1", "UMAPTNK_2")

# Add UMAP reduction to Seurat object
kidney_sle_obj[["umap_tnk"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAPTNK_",
  assay = DefaultAssay(kidney_sle_obj)
)

# Inspect metadata
ncol(kidney_sle_obj) # 58555 cells
table(kidney_sle_obj$cellType) # includes myeloid celltypes - removed in steps below

# Exclude non-T/ILC/NK cell types from object
tnk_celltypes <- unique(as.character(kidney_sle_obj$cellType))
tnk_celltypes <- tnk_celltypes[grepl("^(T|NK)", tnk_celltypes)]
tnk_celltypes
kidney_sle_obj <- subset(kidney_sle_obj, subset = cellType %in% tnk_celltypes)

# Inspect metadata
ncol(kidney_sle_obj) # 34736 cells
table(kidney_sle_obj$cellType)

# Perform typical Seurat workflow of log-normalization, feature selection, and feature scaling
# Reflecting Al Souz et al. methods - "RNA data was normalized by LogNormalize (feature expression measurement divided by the total expression in each cell multiplied by a scale factor of 10,000 then natural log transformed) in Seurat v5. Top 2000 variable genes were determined using FindVariableFeatures with ‘vst’ selection method, and the normalized RNA data was scaled for variable genes using ScaleData."
kidney_sle_obj <- NormalizeData(kidney_sle_obj)
kidney_sle_obj <- FindVariableFeatures(kidney_sle_obj)
kidney_sle_obj <- ScaleData(kidney_sle_obj, features = rownames(kidney_sle_obj))

# Save recreated Seurat object
saveRDS(kidney_sle_obj, '/filepath/fig6/kidney_sle_seurat_obj_recreated.rds')

# Fig. S16C - Annotation dotplot of Tfh/Tph-like cells in Al Souz et al. adult healthy versus SLE kidney sample scRNAseq dataset ----

# Import kidney T/ILC/NK Seurat object from Al Souz et al. (code to recreate object provided above)
kidney_sle_obj <- readRDS('/filepath/fig6/kidney_sle_seurat_obj_recreated.rds')

# Rename gd T cells (to avoid symbol graphical error when saving figure)
kidney_sle_obj$cellType <- gsub(
  "γδ",
  "gd",
  kidney_sle_obj$cellType
)

# Order clusters for dotplot annotation
celltype_order <- unique(kidney_sle_obj$cellType) %>%
  as.character() %>%
  tibble::tibble(cellType = .) %>%
  dplyr::mutate(
    cluster_number = as.numeric(stringr::str_extract(cellType, "\\d+"))
  ) %>%
  dplyr::arrange(cluster_number, cellType) %>%
  dplyr::pull(cellType)

kidney_sle_obj$cellType <- factor(
  kidney_sle_obj$cellType,
  levels = rev(celltype_order)
)

# Assemble annotation dotplot
pdf('/filepath/fig6/fig6_supp/kidney_ln_tnk_annotation_dotplot.pdf',
    height = 8,
    width = 6.75)

DotPlot(
  kidney_sle_obj,
  features = c(
    'CXCR5','PDCD1','ICOS','BTLA','CXCL13','IL21',
    'BCL6','TOX2','POU2AF1','NFATC1','BACH1','GNG4'
  ), 
  group.by = 'cellType',
  dot.scale = 9,
  cluster.idents = FALSE
) +
  scale_color_gradient2(
    low = "#4E79A7",
    mid = "white",
    high = "coral3",
    midpoint = 0
  ) +
  scale_x_discrete(
    labels = function(x) ifelse(x == "GNG4", "***GNG4***", paste0("*", x, "*")) # bold GNG4 gene name
  ) +
  scale_y_discrete(
    labels = function(x) ifelse(x == "T14. Tfh/Tph", "**T14. Tfh/Tph**", x) # bold Tfh/Tph-like cluster name
  ) +
  guides(
    color = guide_colorbar(
      title = "Scaled RNA expression",
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(2.5, "cm"),
      barheight = grid::unit(0.35, "cm")
    ),
    size = guide_legend(
      title = "Percent expressed",
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE,
      keywidth = grid::unit(0.45, "cm"),
      keyheight = grid::unit(0.35, "cm")
    )
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
    axis.text.y = ggtext::element_markdown(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

dev.off()

# Fig. S16D - % GNG4 RNA+ of each subset in Al Souz et al. adult healthy versus SLE kidney sample scRNAseq dataset ----

# Define GNG4 RNA+ cells
DefaultAssay(kidney_sle_obj) <- 'RNA'
kidney_sle_obj$GNG4_pos <- FetchData(kidney_sle_obj, vars = "GNG4")[,1] > 0

# Get metadata
gng4_pct_df <- kidney_sle_obj@meta.data %>%
  dplyr::select(cellType, GNG4_pos) %>%
  dplyr::group_by(cellType) %>%
  dplyr::summarise(
    total_n = n(),
    gng4_pos_n = sum(GNG4_pos),
    pct_gng4_pos = 100 * gng4_pos_n / total_n,
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(pct_gng4_pos))

# Assemble barplot
pdf('/filepath/fig6/fig6_supp/kidney_ln_tnk_pct_gng4pos_by_cluster_barplot.pdf',
    height = 6,
    width = 5.75)
ggplot(gng4_pct_df, aes(x = reorder(cellType, pct_gng4_pos), y = pct_gng4_pos)) +
  geom_bar(stat = "identity", fill = "coral3", color = "black", width = 0.8) +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) ifelse(
      x == "T14. Tfh/Tph",
      "**T14. Tfh/Tph**",
      x
    )
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "" # expression("% " * italic('GNG4') * " RNA+ of cells within cluster")
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.y = ggtext::element_markdown(size = 12),
    axis.text.x = element_text(size = 12)
  )
dev.off()

# Fig. S16E - % GNG4 RNA+ of each subset in Balasubramanian et al. pediatric healthy versus SLE PBMC scRNAseq dataset ----

# Reanalysis of Balasubramanian et al. (2025) Nat Immunol
# Single-cell RNA profiling of blood CD4+ T cells identifies distinct helper and dysfunctional regulatory clusters in children with SLE
# https://pubmed.ncbi.nlm.nih.gov/41120754/
# Retrieved 'sorted.cd4.total.rds' file from Zenodo repository - https://zenodo.org/records/16385794

# Import Seurat object
pbmc_sle_obj <- readRDS('/filepath/fig6/sorted.cd4.total.rds')

# Annotation of activated cTfh-like cluster (based on reanalysis including features in Fig. 6H)
pbmc_sle_obj$clusters_grouped <- gsub(
  "Proliferating",
  "Act. cTfh-like",
  pbmc_sle_obj$clusters_grouped
)

# Define GNG4 RNA positive cells
DefaultAssay(pbmc_sle_obj) <- 'RNA'
pbmc_sle_obj$GNG4_pos <- FetchData(pbmc_sle_obj, vars = "GNG4")[,1] > 0

# Get metadata
gng4_pct_df <- pbmc_sle_obj@meta.data %>%
  dplyr::select(clusters_grouped, GNG4_pos) %>%
  dplyr::group_by(clusters_grouped) %>%
  dplyr::summarise(
    total_n = n(),
    gng4_pos_n = sum(GNG4_pos),
    pct_gng4_pos = 100 * gng4_pos_n / total_n,
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(pct_gng4_pos))

# Assemble barplot
pdf('/filepath/fig6/fig6_supp/pbmc_sle_pct_gng4pos_by_cluster_barplot.pdf',
    height = 6,
    width = 5)
ggplot(gng4_pct_df, aes(x = reorder(clusters_grouped, pct_gng4_pos), y = pct_gng4_pos)) +
  geom_bar(stat = "identity", fill = "coral3", color = "black", width = 0.8) +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) ifelse(
      x == "Act. cTfh-like",
      "**Activated<br>cTfh-like**",
      x
    )
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "" # expression("% " * italic('GNG4') * " RNA+ of cells within cluster")
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text.y = ggtext::element_markdown(size = 12),
    axis.text.x = element_text(size = 12)
  )
dev.off()