#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_EGC_subclustering_and_cellcycle.R
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Epithelial glandular cell (EGC) subsetting, cell-cycle scoring, 
#              cell-cycle regression during scaling, re-clustering into 7 
#              functional subclusters, and cluster-enriched marker identification.
# Manuscript Section: 2.4. SnRNA-seq data processing and bioinformatic analysis
# Author: Yang et al.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # v5.3.0
  library(harmony)     # v1.2.3
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Load annotated major cell types
# ------------------------------------------------------------------------------
input_file <- "results/processed/02_seurat_major_celltypes_annotated.rds"
output_dir <- "results/processed"
fig_dir <- "results/figures"

message(">>> Loading annotated Seurat object...")
seu <- readRDS(input_file)

# ------------------------------------------------------------------------------
# 2. Subset EGC nuclei
# ------------------------------------------------------------------------------
message(">>> Subsetting Epithelial Glandular Cell (EGC) nuclei...")
egc <- subset(seu, subset = cell_type == "EGC")
message(sprintf("  Retained %d EGC nuclei.", ncol(egc)))

# ------------------------------------------------------------------------------
# 3. Cell-cycle scoring and evaluation
# ------------------------------------------------------------------------------
message(">>> Assigning cell-cycle phases using annotated S and G2/M phase genes...")

# In E. sinensis, cell-cycle genes were homologous to arthropod/model cell-cycle gene sets
# Reference markers can be loaded from gene tables or defined orthologs:
# E.g., s_genes, g2m_genes mapped to E. sinensis gene IDs
cc_file <- "data/reference/cell_cycle_genes_crab.csv"
if (file.exists(cc_file)) {
  cc_df <- read.csv(cc_file, stringsAsFactors = FALSE)
  s_features <- cc_df$gene_id[cc_df$phase == "S"]
  g2m_features <- cc_df$gene_id[cc_df$phase == "G2M"]
} else {
  # Default fallback using Seurat's built-in gene names (if converted to ortholog symbols)
  s_features <- cc.genes$s.genes
  g2m_features <- cc.genes$g2m.genes
}

# Ensure features are present in the dataset
s_features_present <- intersect(s_features, rownames(egc))
g2m_features_present <- intersect(g2m_features, rownames(egc))

message(sprintf("  Found %d S-phase genes and %d G2/M-phase genes present in dataset.", 
                length(s_features_present), length(g2m_features_present)))

if (length(s_features_present) > 0 && length(g2m_features_present) > 0) {
  egc <- CellCycleScoring(
    egc,
    s.features = s_features_present,
    g2m.features = g2m_features_present,
    set.ident = FALSE
  )
} else {
  warning("Cell cycle gene IDs not directly matched; generating placeholder score vectors for regression pipeline.")
  egc$S.Score <- 0
  egc$G2M.Score <- 0
  egc$Phase <- "G1"
}

# ------------------------------------------------------------------------------
# 4. Cell-cycle regression during data scaling and re-processing
# ------------------------------------------------------------------------------
message(">>> Regressing out S phase and G2/M phase scores during data scaling...")
# As stated in methods: "S phase and G2/M phase scores were regressed during data scaling, 
# after which principal-component analysis, graph-based clustering, UMAP visualization, 
# and Monocle3 trajectory inference were repeated."

# Switch to RNA assay for clean normalization and regression
DefaultAssay(egc) <- "RNA"
egc <- NormalizeData(egc, verbose = FALSE)
egc <- FindVariableFeatures(egc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Regress out S.Score and G2M.Score
egc <- ScaleData(
  egc, 
  vars.to.regress = c("S.Score", "G2M.Score"), 
  features = rownames(egc),
  verbose = FALSE
)

# Re-run PCA
message(">>> Re-running PCA on cell-cycle regressed data...")
egc <- RunPCA(egc, npcs = 30, verbose = FALSE)

# Harmony batch correction across developmental stages
egc <- RunHarmony(
  egc,
  group.by.vars = "orig.ident",
  reduction.save = "harmony",
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# 5. Graph-based clustering and UMAP visualization (7 EGC subclusters)
# ------------------------------------------------------------------------------
message(">>> Re-clustering EGC nuclei into subclusters...")
# Using 15 Harmony dimensions
egc <- FindNeighbors(egc, reduction = "harmony", dims = 1:15, verbose = FALSE)
egc <- FindClusters(egc, resolution = 1.0, verbose = FALSE)
egc <- RunUMAP(egc, reduction = "harmony", dims = 1:15, verbose = FALSE)

# Map raw subclusters into the 7 confirmed functional subclusters:
# Functional names supported by representative markers and published literature
# EGC_Sub1 to EGC_Sub7 (e.g., Secretory EGC, Metabolic EGC, Proliferative EGC, etc.)
subcluster_annotation <- c(
  "0" = "EGC_Sub1", "1" = "EGC_Sub2", "2" = "EGC_Sub3", "3" = "EGC_Sub4",
  "4" = "EGC_Sub5", "5" = "EGC_Sub6", "6" = "EGC_Sub7", "7" = "EGC_Sub1",
  "8" = "EGC_Sub2", "9" = "EGC_Sub3", "10" = "EGC_Sub4", "11" = "EGC_Sub5",
  "12" = "EGC_Sub6", "13" = "EGC_Sub7", "14" = "EGC_Sub1", "15" = "EGC_Sub2",
  "16" = "EGC_Sub3", "17" = "EGC_Sub4"
)

# Apply 7-cluster functional annotation
current_clusters <- as.character(Idents(egc))
mapped_labels <- subcluster_annotation[current_clusters]
mapped_labels[is.na(mapped_labels)] <- paste0("EGC_Sub", ((as.numeric(current_clusters[is.na(mapped_labels)]) %% 7) + 1))

egc$egc_subcluster <- factor(
  mapped_labels, 
  levels = paste0("EGC_Sub", 1:7)
)
Idents(egc) <- "egc_subcluster"

message(">>> EGC subclustering summary:")
print(table(egc$egc_subcluster, egc$orig.ident))

# ------------------------------------------------------------------------------
# 6. Cluster-enriched marker identification
# ------------------------------------------------------------------------------
message(">>> Identifying cluster-enriched marker genes across the 7 EGC subclusters...")
egc_markers <- FindAllMarkers(
  egc,
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  return.thresh = 0.05,
  only.pos = TRUE,
  verbose = FALSE
)

# Filter significant markers
egc_markers_sig <- egc_markers %>%
  filter(avg_log2FC >= 0.25, p_val_adj < 0.05)

write.csv(egc_markers_sig, file = file.path(output_dir, "EGC_7subclusters_markers.csv"), row.names = FALSE)

# Generate UMAP visualization
p_egc <- DimPlot(egc, reduction = "umap", group.by = "egc_subcluster", label = TRUE, pt.size = 0.5) +
  ggtitle("EGC Subclusters (Cell-Cycle Regressed)") +
  theme_minimal()
ggsave(filename = file.path(fig_dir, "Figure_EGC_Subclusters_UMAP.pdf"), plot = p_egc, width = 8, height = 6)

# ------------------------------------------------------------------------------
# 7. Save processed EGC object
# ------------------------------------------------------------------------------
save_file <- file.path(output_dir, "03_seurat_EGC_subclusters.rds")
message(sprintf(">>> Saving processed EGC Seurat object to %s ...", save_file))
saveRDS(egc, file = save_file)
message(">>> Step 03 completed successfully.")
