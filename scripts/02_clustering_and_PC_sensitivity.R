#!/usr/bin/env Rscript
# ==============================================================================
# Script: 02_clustering_and_PC_sensitivity.R
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Graph-based clustering with 10 Harmony PCs (res = 0.5), 
#              principal component (PC) sensitivity analysis (15, 20, 30, 50 PCs; 
#              Figure S1), candidate marker identification via bimod test, and
#              major cell-type annotation (EGC, SSC, IC, PC).
# Manuscript Section: 2.4. SnRNA-seq data processing and bioinformatic analysis
# Author: Yang et al.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # v5.3.0
  library(harmony)     # v1.2.3
  library(mclust)      # for adjustedRandIndex (ARI)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. Load integrated data
# ------------------------------------------------------------------------------
input_file <- "results/processed/01_seurat_harmony_integrated.rds"
output_dir <- "results/processed"
fig_dir <- "results/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

message(">>> Loading Harmony-integrated Seurat object...")
seu <- readRDS(input_file)

# ------------------------------------------------------------------------------
# 2. Reference clustering using the first 10 Harmony-corrected dimensions
# ------------------------------------------------------------------------------
message(">>> Performing reference clustering using 10 Harmony PCs (resolution = 0.5)...")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:10, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:10, verbose = FALSE)

# Store baseline cluster assignments
seu$ref_clusters_10pc <- Idents(seu)

# ------------------------------------------------------------------------------
# 3. Major cell type marker detection via bimod test and annotation
# ------------------------------------------------------------------------------
message(">>> Identifying cluster markers using the bimod likelihood-ratio test...")
# Candidate markers: expressed in >10% nuclei, avg_log2FC >= 0.25, p_val_adj < 0.05
markers_all <- FindAllMarkers(
  seu,
  assay = "SCT",
  test.use = "bimod",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  return.thresh = 0.05,
  only.pos = TRUE,
  verbose = FALSE
)

# Filter candidate markers
candidate_markers <- markers_all %>%
  filter(avg_log2FC >= 0.25, p_val_adj < 0.05)

write.csv(candidate_markers, file = file.path(output_dir, "major_cluster_candidate_markers_bimod.csv"), row.names = FALSE)

# Cluster mapping to major cell types based on canonical marker expression and histology:
# EGC: Epithelial glandular cells
# SSC: Stromal support cells
# IC:  Immune cells
# PC:  Putative perisinusoidal cells
# (Note: Map cluster IDs to cell types according to reference marker expression profiles)
cluster_to_celltype <- c(
  "0" = "EGC",
  "1" = "EGC",
  "2" = "EGC",
  "3" = "SSC",
  "4" = "EGC",
  "5" = "IC",
  "6" = "PC"
)

# Apply annotations
seu$cell_type <- unname(cluster_to_celltype[as.character(Idents(seu))])
seu$cell_type <- factor(seu$cell_type, levels = c("EGC", "SSC", "IC", "PC"))
Idents(seu) <- "cell_type"
seu$ref_major_celltype <- seu$cell_type

# ------------------------------------------------------------------------------
# 4. Principal Component (PC) Sensitivity Analysis (Figure S1)
# ------------------------------------------------------------------------------
message(">>> Conducting PC sensitivity analysis across 15, 20, 30, and 50 Harmony dimensions...")
# Evaluate sensitivity against original 10-PC reference using ARI and % concordance
test_pcs <- c(15, 20, 30, 50)
sensitivity_metrics <- data.frame(
  PC_dimensions = c(10, test_pcs),
  Adjusted_Rand_Index = c(1.0, rep(NA, length(test_pcs))),
  Concordance_Percent = c(100.0, rep(NA, length(test_pcs)))
)

sensitivity_umap_list <- list()
sensitivity_umap_list[["PC10"]] <- DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE) +
  ggtitle("Reference: 10 PCs") + theme_minimal()

for (i in seq_along(test_pcs)) {
  npc <- test_pcs[i]
  message(sprintf("  Evaluating %d Harmony PCs...", npc))
  
  # Recalculate neighbors, clusters, and UMAP
  temp_seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npc, verbose = FALSE)
  temp_seu <- FindClusters(temp_seu, resolution = 0.5, verbose = FALSE)
  temp_seu <- RunUMAP(temp_seu, reduction = "harmony", dims = 1:npc, verbose = FALSE)
  
  new_clusters <- temp_seu$seurat_clusters
  ref_labels <- seu$ref_major_celltype
  
  # Compute Adjusted Rand Index (ARI) comparing clustering partitions
  ari_val <- adjustedRandIndex(ref_labels, new_clusters)
  
  # Assign each new cluster the most frequent original major-cell-type label among its nuclei
  cluster_mapping <- table(new_clusters, ref_labels)
  cluster_assigned_label <- apply(cluster_mapping, 1, function(row) names(which.max(row)))
  
  projected_labels <- cluster_assigned_label[as.character(new_clusters)]
  concordance_rate <- mean(projected_labels == ref_labels) * 100
  
  sensitivity_metrics[i + 1, "Adjusted_Rand_Index"] <- round(ari_val, 4)
  sensitivity_metrics[i + 1, "Concordance_Percent"] <- round(concordance_rate, 2)
  
  temp_seu$projected_celltype <- factor(projected_labels, levels = levels(ref_labels))
  sensitivity_umap_list[[paste0("PC", npc)]] <- DimPlot(
    temp_seu, 
    reduction = "umap", 
    group.by = "projected_celltype", 
    label = TRUE
  ) + ggtitle(sprintf("%d PCs (Concordance: %.1f%%)", npc, concordance_rate)) + theme_minimal()
}

write.csv(sensitivity_metrics, file = file.path(output_dir, "PC_sensitivity_metrics_FigS1.csv"), row.names = FALSE)
message(">>> PC Sensitivity Metrics:")
print(sensitivity_metrics)

# Assemble Figure S1
p_s1 <- (sensitivity_umap_list[["PC10"]] + sensitivity_umap_list[["PC15"]]) /
        (sensitivity_umap_list[["PC20"]] + sensitivity_umap_list[["PC30"]]) /
        (sensitivity_umap_list[["PC50"]] + plot_spacer()) +
        plot_annotation(title = "Figure S1: Sensitivity Analysis across Harmony PC Dimensions")

ggsave(filename = file.path(fig_dir, "Figure_S1_PC_Sensitivity.pdf"), plot = p_s1, width = 12, height = 15)
ggsave(filename = file.path(fig_dir, "Figure_S1_PC_Sensitivity.png"), plot = p_s1, width = 12, height = 15, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Save annotated Seurat object
# ------------------------------------------------------------------------------
save_file <- file.path(output_dir, "02_seurat_major_celltypes_annotated.rds")
message(sprintf(">>> Saving annotated Seurat object to %s ...", save_file))
saveRDS(seu, file = save_file)
message(">>> Step 02 completed successfully.")
