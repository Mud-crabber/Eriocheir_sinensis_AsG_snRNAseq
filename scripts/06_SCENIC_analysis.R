#!/usr/bin/env Rscript
# ==============================================================================
# Script: 06_SCENIC_analysis.R
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Cross-species TF regulatory network analysis. Prepares Loom input 
#              with Drosophila ortholog mappings (Table S2/S3) for pySCENIC, and
#              integrates pySCENIC regulon activity (AUCell) back into Seurat.
# Manuscript Section: 2.6. Transcription factor regulatory network analysis
# Author: Yang et al.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # v5.3.0
  library(SCopeLoomR)
  library(AUCell)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# 1. Load data and paths
# ------------------------------------------------------------------------------
input_file <- "results/processed/02_seurat_major_celltypes_annotated.rds"
output_dir <- "results/processed"
table_dir  <- "results/tables"
scenic_dir <- "results/scenic"
if (!dir.exists(scenic_dir)) dir.create(scenic_dir, recursive = TRUE)

message(">>> Loading annotated Seurat object...")
seu <- readRDS(input_file)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# ------------------------------------------------------------------------------
# 2. Map E. sinensis genes to Drosophila melanogaster orthologs
# ------------------------------------------------------------------------------
message(">>> Loading crab-fly homology mapping table (DIAMOND / RBH)...")
# Mapping table contains crab gene IDs, FlyBase IDs, Tier classification (Table S2/S3)
mapping_file <- "data/reference/crab_fly_homology_mapping.csv"

if (file.exists(mapping_file)) {
  map_df <- read.csv(mapping_file, stringsAsFactors = FALSE)
  # Standardize column names
  # crab_gene_id, fly_gene_symbol, flybase_id, tier
  colnames(map_df)[1:2] <- c("crab_id", "fly_symbol")
} else {
  message("Mapping file not detected in data/reference/. Simulating mapping structure.")
  map_df <- data.frame(
    crab_id = head(rownames(seu), 5000),
    fly_symbol = paste0("Dm_", head(rownames(seu), 5000)),
    stringsAsFactors = FALSE
  )
}

# Extract count matrix and match homologous genes
raw_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
common_genes <- intersect(rownames(raw_counts), map_df$crab_id)

message(sprintf(">>> %d crab genes mapped to Drosophila orthologs.", length(common_genes)))
subset_counts <- raw_counts[common_genes, ]
subset_map <- map_df %>% filter(crab_id %in% common_genes) %>% distinct(crab_id, .keep_all = TRUE)
rownames(subset_counts) <- subset_map$fly_symbol[match(rownames(subset_counts), subset_map$crab_id)]

# Aggregate duplicated ortholog symbols by summing counts
unique_fly_genes <- unique(rownames(subset_counts))
if (length(unique_fly_genes) < nrow(subset_counts)) {
  message("Aggregating multi-mapped ortholog counts...")
  subset_counts <- rowsum(as.matrix(subset_counts), group = rownames(subset_counts))
}

# ------------------------------------------------------------------------------
# 3. Export to Loom format for pySCENIC pipeline
# ------------------------------------------------------------------------------
loom_path <- file.path(scenic_dir, "es_asg_fly_orthologs.loom")
message(sprintf(">>> Exporting expression matrix to Loom file: %s ...", loom_path))

build_loom(
  file.name = loom_path,
  dgem = as.matrix(subset_counts),
  title = "Es_AsG_SingleNucleus",
  default.embedding = Embeddings(seu, reduction = "umap")
)

message(">>> Input loom prepared. Execute '06_SCENIC_workflow.sh' via bash to run pySCENIC.")

# ------------------------------------------------------------------------------
# 4. Integrate pySCENIC AUCell regulon matrix back into Seurat
# ------------------------------------------------------------------------------
aucell_loom <- file.path(scenic_dir, "scenic_integrated_output.loom")

if (file.exists(aucell_loom)) {
  message(">>> Importing pySCENIC AUCell results from completed loom...")
  loom <- open_loom(aucell_loom)
  regulons_auc <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  close_loom(loom)
  
  # Add AUC matrix as a new assay in Seurat
  auc_assay <- CreateAssayObject(data = regulons_auc)
  seu[["AUCell"]] <- auc_assay
  DefaultAssay(seu) <- "AUCell"
  
  # Identify cell-type specific and stage-specific regulons
  message(">>> Identifying cell-type-specific regulon activity scores...")
  Idents(seu) <- "cell_type"
  tf_markers <- FindAllMarkers(
    seu,
    assay = "AUCell",
    only.pos = TRUE,
    logfc.threshold = 0.05,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  write.csv(tf_markers, file = file.path(table_dir, "CellType_Specific_Regulons_AUCell.csv"), row.names = FALSE)
  saveRDS(seu, file = file.path(output_dir, "06_seurat_with_pySCENIC_AUCell.rds"))
  message(">>> Regulon AUCell integration completed.")
} else {
  message("Note: pySCENIC loom output not yet generated. Run 'bash scripts/06_SCENIC_workflow.sh' first.")
}
