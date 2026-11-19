#!/usr/bin/env bash
# ==============================================================================
# Script: 06_SCENIC_workflow.sh
# Project: Single-nucleus transcriptomics of E. sinensis accessory sex gland
# Description: Cross-species pySCENIC (v0.12.1) command-line execution pipeline.
#              Runs GRNBoost2, RcisTarget motif enrichment, and AUCell scoring.
# Manuscript Section: 2.6. Transcription factor regulatory network analysis
# Author: Yang et al.
# ==============================================================================

set -e

# Working directories
DATA_DIR="results/scenic"
REF_DIR="data/reference/scenic_dmel"
THREADS=16

# Reference database files for Drosophila melanogaster (FlyBase / Aerts lab)
# Downloaded from https://resources.aertslab.org/cistarget/
TF_LIST="${REF_DIR}/allTFs_dmel.txt"
MOTIF_RANKINGS="${REF_DIR}/dm6_v10_refseq-r6.01_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
MOTIF_ANNOTATION="${REF_DIR}/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"

INPUT_LOOM="${DATA_DIR}/es_asg_fly_orthologs.loom"
ADJACENCIES="${DATA_DIR}/adjacencies.tsv"
REGULONS="${DATA_DIR}/regulons.csv"
OUTPUT_LOOM="${DATA_DIR}/scenic_integrated_output.loom"

echo "=== [Step 1/3] Running GRNBoost2 (Inference of co-expression modules) ==="
pyscenic grn \
    --num_workers ${THREADS} \
    --output ${ADJACENCIES} \
    --method grnboost2 \
    ${INPUT_LOOM} \
    ${TF_LIST}

echo "=== [Step 2/3] Running RcisTarget (Motif enrichment & regulon pruning) ==="
pyscenic ctx \
    ${ADJACENCIES} \
    ${MOTIF_RANKINGS} \
    --annotations_fname ${MOTIF_ANNOTATION} \
    --expression_mtx_fname ${INPUT_LOOM} \
    --output ${REGULONS} \
    --num_workers ${THREADS}

echo "=== [Step 3/3] Running AUCell (Quantifying regulon activity per nucleus) ==="
pyscenic aucell \
    ${INPUT_LOOM} \
    ${REGULONS} \
    --output ${OUTPUT_LOOM} \
    --num_workers ${THREADS}

echo "=== pySCENIC pipeline completed successfully. Output saved to: ${OUTPUT_LOOM} ==="
