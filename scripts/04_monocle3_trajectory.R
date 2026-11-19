#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_monocle3_trajectory.R
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Single-nucleus pseudotime trajectory inference using Monocle3 (v1.4.27)
#              with cell-cycle regressed EGC nuclei. Evaluates trajectory dynamics
#              across AsG developmental stages (orig.ident: ASG1, ASG2, ASG3).
# Manuscript Section: 2.5. Pseudotime trajectory analysis
# Author: Yang et al.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # v5.3.0
  library(monocle3)    # v1.4.27
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. Load cell-cycle regressed EGC Seurat object
# ------------------------------------------------------------------------------
input_file <- "results/processed/03_seurat_EGC_subclusters.rds"
output_dir <- "results/processed"
fig_dir <- "results/figures"

message(">>> Loading cell-cycle regressed EGC Seurat object...")
egc <- readRDS(input_file)

# In Seurat v5, merge expression layers across samples using JoinLayers
message(">>> Merging expression layers across samples using JoinLayers (Seurat v5)...")
DefaultAssay(egc) <- "RNA"
egc <- JoinLayers(egc)

# ------------------------------------------------------------------------------
# 2. Extract counts and construct Monocle3 cell_data_set (CDS)
# ------------------------------------------------------------------------------
message(">>> Extracting raw count matrix and metadata for Monocle3...")
# Extract raw counts
raw_counts <- GetAssayData(egc, assay = "RNA", layer = "counts")

# Nuclei metadata (cell_metadata)
cell_metadata <- egc@meta.data

# Gene metadata (gene_metadata)
gene_metadata <- data.frame(
  gene_short_name = rownames(raw_counts),
  row.names = rownames(raw_counts)
)

# Construct Monocle3 CDS object
message(">>> Constructing Monocle3 cell_data_set using new_cell_data_set()...")
cds <- new_cell_data_set(
  expression_data = raw_counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# ------------------------------------------------------------------------------
# 3. Preprocess and infer developmental trajectories
# ------------------------------------------------------------------------------
message(">>> Preprocessing CDS using preprocess_cds(num_dim = 50)...")
cds <- preprocess_cds(cds, num_dim = 50)

# Transfer Seurat UMAP embeddings to Monocle3 CDS for visual and topological consistency
message(">>> Transferring Seurat UMAP embeddings into Monocle3 CDS...")
reducedDims(cds)$UMAP <- Embeddings(egc, reduction = "umap")

# Invert / learn trajectory graph without partitioning as stated in manuscript:
# "Cell developmental trajectories were inferred using learn_graph (use_partition = FALSE)"
message(">>> Learning trajectory graph (use_partition = FALSE)...")
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# ------------------------------------------------------------------------------
# 4. Order cells and calculate pseudotime
# ------------------------------------------------------------------------------
message(">>> Ordering cells along pseudotime...")
# Root cells are defined based on the pre-developmental stage (ASG1)
# or early-stage progenitor subcluster
asg1_cells <- which(colData(cds)$orig.ident == "ASG1")
root_node <- names(which.max(table(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[asg1_cells])))

if (!is.null(root_node)) {
  cds <- order_cells(cds, root_pr_nodes = root_node)
} else {
  # Fallback to interactive or default root selection
  cds <- order_cells(cds)
}

# Add pseudotime values back to Seurat metadata
message(">>> Adding pseudotime values to Seurat metadata for downstream visualization...")
egc$pseudotime <- pseudotime(cds)

# ------------------------------------------------------------------------------
# 5. Spatial autocorrelation / trajectory-dependent gene discovery (graph_test)
# ------------------------------------------------------------------------------
message(">>> Identifying trajectory-dependent genes using graph_test (Moran's I)...")
# Run graph_test on variable genes
var_genes <- VariableFeatures(egc)
if (length(var_genes) == 0) var_genes <- head(rownames(raw_counts), 2000)

pr_graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
pr_graph_sig <- pr_graph_test_res %>%
  filter(status == "OK", q_value < 0.05) %>%
  arrange(desc(morans_I))

write.csv(pr_graph_sig, file = file.path(output_dir, "Monocle3_trajectory_genes_moransI.csv"), row.names = TRUE)

# ------------------------------------------------------------------------------
# 6. Visualization across developmental stages (orig.ident: ASG1, ASG2, ASG3)
# ------------------------------------------------------------------------------
message(">>> Generating pseudotime visualizations...")

# 1) Pseudotime FeaturePlot
p_pseudo <- FeaturePlot(egc, features = "pseudotime", cols = c("#440154FF", "#21908CFF", "#FDE725FF")) +
  ggtitle("EGC Nuclei Pseudotime Trajectory") + theme_minimal()

# 2) Pseudotime comparisons among developmental stages (orig.ident)
p_box <- ggplot(egc@meta.data, aes(x = orig.ident, y = pseudotime, fill = orig.ident)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Pseudotime Distribution Across Stages",
    x = "Developmental Stage (orig.ident)",
    y = "Pseudotime"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine plots
p_combined <- p_pseudo + p_box
ggsave(filename = file.path(fig_dir, "Figure_Pseudotime_Trajectory.pdf"), plot = p_combined, width = 12, height = 5)
ggsave(filename = file.path(fig_dir, "Figure_Pseudotime_Trajectory.png"), plot = p_combined, width = 12, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 7. Save updated Seurat object and Monocle3 CDS
# ------------------------------------------------------------------------------
saveRDS(cds, file = file.path(output_dir, "04_monocle3_cds_trajectory.rds"))
saveRDS(egc, file = file.path(output_dir, "04_seurat_EGC_with_pseudotime.rds"))
message(">>> Step 04 completed successfully.")
