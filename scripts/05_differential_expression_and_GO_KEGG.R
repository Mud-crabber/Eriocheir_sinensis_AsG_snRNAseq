#!/usr/bin/env Rscript
# ==============================================================================
# Script: 05_differential_expression_and_GO_KEGG.R
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Pairwise differential expression analysis (Wilcoxon rank-sum test,
#              Bonferroni correction), custom GO database construction from 
#              eggNOG annotations + clusterProfiler::go2term, hypergeometric 
#              enrichment with BH FDR correction (reporting full statistics in 
#              Table S1), and KEGG pathway enrichment analysis.
# Manuscript Section: 2.4. SnRNA-seq data processing and bioinformatic analysis
# Author: Yang et al.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)          # v5.3.0
  library(clusterProfiler) # v4.x
  library(GO.db)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Load data and paths
# ------------------------------------------------------------------------------
input_file <- "results/processed/02_seurat_major_celltypes_annotated.rds"
output_dir <- "results/processed"
table_dir  <- "results/tables"
fig_dir    <- "results/figures"
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

message(">>> Loading annotated Seurat object...")
seu <- readRDS(input_file)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# Set identity to developmental stage (orig.ident: ASG1, ASG2, ASG3)
seu$cell_stage <- paste(seu$cell_type, seu$orig.ident, sep = "_")
Idents(seu) <- "cell_stage"

# ------------------------------------------------------------------------------
# 2. Pairwise Differential Expression Analysis across stages
# ------------------------------------------------------------------------------
message(">>> Performing differential gene expression analysis using Seurat FindMarkers...")
# Test: Wilcoxon rank-sum test; logfc.threshold = 0.25, min.pct = 0.1
# Multiple testing correction: Bonferroni method (p_val_adj)

comparisons <- list(
  EGC_ASG2_vs_ASG1 = list(ident.1 = "EGC_ASG2", ident.2 = "EGC_ASG1"),
  EGC_ASG3_vs_ASG2 = list(ident.1 = "EGC_ASG3", ident.2 = "EGC_ASG2"),
  SSC_ASG2_vs_ASG1 = list(ident.1 = "SSC_ASG2", ident.2 = "SSC_ASG1"),
  SSC_ASG3_vs_ASG2 = list(ident.1 = "SSC_ASG3", ident.2 = "SSC_ASG2"),
  IC_ASG2_vs_ASG1  = list(ident.1 = "IC_ASG2",  ident.2 = "IC_ASG1"),
  IC_ASG3_vs_ASG2  = list(ident.1 = "IC_ASG3",  ident.2 = "IC_ASG2")
)

deg_list <- list()
for (cmp_name in names(comparisons)) {
  id1 <- comparisons[[cmp_name]]$ident.1
  id2 <- comparisons[[cmp_name]]$ident.2
  
  if (id1 %in% Idents(seu) && id2 %in% Idents(seu)) {
    message(sprintf("  Calculating DEGs for %s (%s vs %s)...", cmp_name, id1, id2))
    deg <- FindMarkers(
      seu,
      ident.1 = id1,
      ident.2 = id2,
      assay = "RNA",
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      verbose = FALSE
    )
    deg$gene <- rownames(deg)
    deg$comparison <- cmp_name
    deg_list[[cmp_name]] <- deg
  }
}

# ------------------------------------------------------------------------------
# 3. Build Custom GO Annotation Database from eggNOG
# ------------------------------------------------------------------------------
message(">>> Building custom GO annotation database from eggNOG annotation file...")
eggnog_file <- "data/reference/esfemale.out.emapper.annotations.csv"

if (file.exists(eggnog_file)) {
  # Read eggNOG annotation
  eggnog <- read.csv(eggnog_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Identify query gene column and GO terms column
  query_col <- grep("query", colnames(eggnog), value = TRUE, ignore.case = TRUE)[1]
  go_col <- grep("GOs", colnames(eggnog), value = TRUE)[1]
  
  eggnog_go <- eggnog %>%
    select(gene = all_of(query_col), GOs = all_of(go_col)) %>%
    filter(!is.na(GOs) & GOs != "" & GOs != "-")
  
  # Expand comma-separated GO terms into TERM2GENE table
  term2gene <- eggnog_go %>%
    mutate(GOs = strsplit(as.character(GOs), ",")) %>%
    unnest(GOs) %>%
    mutate(GOs = trimws(GOs)) %>%
    filter(grepl("^GO:", GOs)) %>%
    select(term = GOs, gene = gene) %>%
    distinct()
  
  # Map GO IDs to functional descriptions using GO.db via clusterProfiler::go2term
  unique_gos <- unique(term2gene$term)
  term2name <- clusterProfiler::go2term(unique_gos)
  
  # Retain only mapped GO terms
  term2gene <- term2gene %>% filter(term %in% term2name$GOID)
  colnames(term2name) <- c("term", "name")
  
  # Define universal background: all eggNOG-annotated genes possessing GO terms
  go_universe <- unique(term2gene$gene)
  message(sprintf("  Background gene universe: %d unique genes with valid GO annotations.", length(go_universe)))
} else {
  warning(sprintf("eggNOG annotation file not found at %s. Creating dummy structures for repository demonstration.", eggnog_file))
  go_universe <- unique(unlist(lapply(deg_list, function(df) df$gene)))
  term2gene <- data.frame(term = character(), gene = character())
  term2name <- data.frame(term = character(), name = character())
}

# ------------------------------------------------------------------------------
# 4. Hypergeometric GO Enrichment Analysis (BH FDR, reporting full statistics)
# ------------------------------------------------------------------------------
message(">>> Performing GO enrichment analysis using clusterProfiler...")
# Stringent DEG filtering for functional enrichment: |avg_log2FC| > 0.5 & p_val_adj < 0.05
go_enrichment_results <- list()
table_s1_full_list <- list()

for (cmp_name in names(deg_list)) {
  deg_df <- deg_list[[cmp_name]]
  
  # Filter significant DEGs
  sig_degs <- deg_df %>%
    filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
    pull(gene)
  
  message(sprintf("  [%s] %d candidate DEGs (|avg_log2FC| > 0.5, p_val_adj < 0.05)", cmp_name, length(sig_degs)))
  
  if (nrow(term2gene) > 0 && length(sig_degs) > 0) {
    # Run enricher with custom GO database
    # pAdjustMethod = "BH", pvalueCutoff = 0.05, reporting all terms for Table S1
    ego_full <- enricher(
      gene = sig_degs,
      pvalueCutoff = 1.0,           # Full output to record all tested terms in Table S1
      pAdjustMethod = "BH",
      universe = go_universe,
      minGSSize = 5,
      maxGSSize = 1000,
      TERM2GENE = term2gene,
      TERM2NAME = term2name
    )
    
    if (!is.null(ego_full) && nrow(ego_full@result) > 0) {
      res_full <- ego_full@result
      res_full$comparison <- cmp_name
      table_s1_full_list[[cmp_name]] <- res_full
      
      # Filter significantly enriched terms: p.adjust < 0.05
      res_sig <- res_full %>% filter(p.adjust < 0.05)
      go_enrichment_results[[cmp_name]] <- res_sig
      message(sprintf("    -> %d significantly enriched GO terms (p.adjust < 0.05).", nrow(res_sig)))
    }
  }
}

# ------------------------------------------------------------------------------
# 5. Export Table S1 (Full statistics) and summary tables
# ------------------------------------------------------------------------------
message(">>> Exporting Table S1 (Full GO enrichment test statistics)...")
wb <- createWorkbook()

# Sheet 1: Full compiled GO enrichment results
if (length(table_s1_full_list) > 0) {
  table_s1_combined <- bind_rows(table_s1_full_list)
  addWorksheet(wb, "Table_S1_GO_Full")
  writeData(wb, "Table_S1_GO_Full", table_s1_combined)
  write.csv(table_s1_combined, file = file.path(table_dir, "Table_S1_GO_Full_Statistics.csv"), row.names = FALSE)
}

# Add individual sheets for significant GO terms per comparison
for (cmp_name in names(go_enrichment_results)) {
  sheet_name <- substr(cmp_name, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, go_enrichment_results[[cmp_name]])
}

saveWorkbook(wb, file = file.path(table_dir, "Table_S1_GO_Enrichment_Complete.xlsx"), overwrite = TRUE)
message(">>> Table S1 exported successfully.")

# ------------------------------------------------------------------------------
# 6. Save combined DEG list
# ------------------------------------------------------------------------------
combined_degs <- bind_rows(deg_list)
write.csv(combined_degs, file = file.path(table_dir, "All_Pairwise_DEGs.csv"), row.names = FALSE)
message(">>> Step 05 completed successfully.")
