#!/usr/bin/env Rscript
# ==============================================================================
# Script: 01_QC_and_harmony_integration.R
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Quality control, doublet removal (DoubletFinder v2.0.6), 
#              SCTransform normalization, and batch correction via Harmony (v1.2.3).
# Manuscript Section: 2.4. SnRNA-seq data processing and bioinformatic analysis
# Author: Yang et al.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)          # v5.3.0
  library(DoubletFinder)   # v2.0.6
  library(harmony)         # v1.2.3
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Define paths and sample metadata
# ------------------------------------------------------------------------------
# Set input and output directories (relative paths for reproducibility)
data_dir <- "data/raw"
output_dir <- "results/processed"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Samples representing three developmental stages:
# ASG1: Pre-developmental stage
# ASG2: Early-developmental stage
# ASG3: Rapid-developmental stage
sample_info <- list(
  ASG1 = file.path(data_dir, "sample_ASG1/filtered_feature_bc_matrix"),
  ASG2 = file.path(data_dir, "sample_ASG2/filtered_feature_bc_matrix"),
  ASG3 = file.path(data_dir, "sample_ASG3/filtered_feature_bc_matrix")
)

# ------------------------------------------------------------------------------
# 2. Load data, QC filtering, and DoubletFinder per sample
# ------------------------------------------------------------------------------
processed_seurat_list <- list()

for (sample_id in names(sample_info)) {
  message(sprintf(">>> Processing sample: %s ...", sample_id))
  matrix_path <- sample_info[[sample_id]]
  
  if (!dir.exists(matrix_path)) {
    warning(sprintf("Directory not found: %s. Please verify data path.", matrix_path))
    next
  }
  
  # Load 10x Cell Ranger output
  counts <- Read10X(data.dir = matrix_path)
  
  # Create Seurat object with min.cells = 3 and min.features = 200
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  # Standard single-nucleus QC filtering as described in manuscript:
  # nCount_RNA < 10,000 and 200 < nFeature_RNA < 3,000
  n_before <- ncol(seu)
  seu <- subset(
    seu, 
    subset = nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 3000
  )
  n_after_qc <- ncol(seu)
  message(sprintf("  [%s] QC Filtering: %d -> %d nuclei retained.", sample_id, n_before, n_after_qc))
  
  # ----------------------------------------------------------------------------
  # Doublet identification and exclusion using DoubletFinder (v2.0.6)
  # ----------------------------------------------------------------------------
  message(sprintf("  [%s] Running DoubletFinder...", sample_id))
  
  # Standard pre-processing required for DoubletFinder
  seu_df <- NormalizeData(seu, verbose = FALSE)
  seu_df <- FindVariableFeatures(seu_df, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu_df <- ScaleData(seu_df, verbose = FALSE)
  seu_df <- RunPCA(seu_df, npcs = 30, verbose = FALSE)
  
  # Parameter sweep to determine optimal pK
  sweep_res <- paramSweep_v3(seu_df, PCs = 1:20, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$MeanBC)]))
  message(sprintf("  [%s] Optimal pK identified: %f", sample_id, optimal_pK))
  
  # Estimate homotypic doublet proportion (assuming ~5% doublet formation rate)
  annotations <- seu_df$seurat_clusters
  homotypic_prop <- if (!is.null(annotations)) modelHomotypic(annotations) else 0.05
  nExp_poi <- round(0.05 * ncol(seu_df))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  
  # Run DoubletFinder
  seu_df <- doubletFinder_v3(
    seu_df,
    PCs = 1:20,
    pN = 0.25,
    pK = optimal_pK,
    nExp = nExp_poi_adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  
  # Identify doublet classification column
  df_class_col <- grep("^DF.classifications", colnames(seu_df@meta.data), value = TRUE)[1]
  seu$DF_classification <- seu_df@meta.data[[df_class_col]]
  
  # Exclude predicted doublets
  seu_singlet <- subset(seu, subset = DF_classification == "Singlet")
  n_singlets <- ncol(seu_singlet)
  message(sprintf("  [%s] Doublet exclusion: %d doublets removed, %d singlets retained.", 
                  sample_id, n_after_qc - n_singlets, n_singlets))
  
  processed_seurat_list[[sample_id]] <- seu_singlet
}

# ------------------------------------------------------------------------------
# 3. Merge samples and batch correction via Harmony (v1.2.3)
# ------------------------------------------------------------------------------
message(">>> Merging samples into a combined Seurat object...")
combined_seurat <- merge(
  x = processed_seurat_list[[1]],
  y = processed_seurat_list[-1],
  add.cell.ids = names(processed_seurat_list),
  project = "Es_AsG"
)

# Set stage factor levels
combined_seurat$stage <- factor(
  combined_seurat$orig.ident,
  levels = c("ASG1", "ASG2", "ASG3"),
  labels = c("Pre-developmental", "Early-developmental", "Rapid-developmental")
)

message(">>> Running SCTransform normalization...")
# Counts were normalized using SCTransform in Seurat
combined_seurat <- SCTransform(
  combined_seurat,
  assay = "RNA",
  variable.features.n = 3000,
  verbose = FALSE
)

message(">>> Running PCA (npcs = 50)...")
combined_seurat <- RunPCA(combined_seurat, assay = "SCT", npcs = 50, verbose = FALSE)

message(">>> Performing batch correction using Harmony (v1.2.3)...")
# Harmony batch correction across orig.ident (sample/batch)
combined_seurat <- RunHarmony(
  combined_seurat,
  group.by.vars = "orig.ident",
  assay.use = "SCT",
  reduction.save = "harmony",
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# 4. Save intermediate integrated object
# ------------------------------------------------------------------------------
save_file <- file.path(output_dir, "01_seurat_harmony_integrated.rds")
message(sprintf(">>> Saving integrated Seurat object to %s ...", save_file))
saveRDS(combined_seurat, file = save_file)
message(">>> Step 01 completed successfully.")
