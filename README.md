# Single-Nucleus RNA Sequencing (snRNA-seq) Analysis of Chinese Mitten Crab (*Eriocheir sinensis*) Accessory Sex Gland

This repository provides the complete, modular, and reproducible bioinformatic analysis workflow for the snRNA-seq study of the Chinese mitten crab (*Eriocheir sinensis*) accessory sex gland (AsG) across different developmental stages (pre-developmental, early-developmental, and rapid-developmental stages).

The analytical methodology implemented here strictly adheres to the procedures described in Section **2.4–2.7** of the manuscript.

---

## Overview of the Workflow

```
Raw FASTQ (10x Genomics)
        │
        ▼ Cell Ranger (v9.0.1; STAR alignment against E. sinensis GCA_013436485.1)
Gene-Barcode Count Matrices
        │
        ├─► 01_QC_and_harmony_integration.R
        │     - Filtering: min.cells = 3, min.features = 200, nCount < 10,000, 200 < nFeature < 3,000
        │     - Doublet exclusion: DoubletFinder (v2.0.6)
        │     - Normalization & Batch Correction: SCTransform, PCA (npcs = 50), Harmony (v1.2.3)
        │
        ├─► 02_clustering_and_PC_sensitivity.R
        │     - Graph clustering: 10 Harmony PCs, resolution = 0.5 (UMAP visualization)
        │     - PC Sensitivity Analysis: evaluated across 15, 20, 30, and 50 PCs (ARI & % concordance, Figure S1)
        │     - Major Cell Types: EGC, SSC, IC, PC identified via FindAllMarkers (bimod test)
        │
        ├─► 03_EGC_subclustering_and_cellcycle.R
        │     - Subclustering of EGC nuclei into 7 functional subclusters (dims 1:15)
        │     - Cell-cycle phase scoring (S and G2/M phases) and regression during data scaling
        │     - Assessment of cell cycle influence on clustering and marker programs
        │
        ├─► 04_monocle3_trajectory.R
        │     - Monocle3 (v1.4.27) trajectory inference: raw counts extraction, preprocess_cds (num_dim = 50)
        │     - Trajectory graph learning (learn_graph, use_partition = FALSE) and pseudotime ordering
        │     - Cell cycle-regressed trajectory validation
        │
        ├─► 05_differential_expression_and_GO_KEGG.R
        │     - Pairwise DEG identification across developmental stages (FindMarkers; Wilcoxon rank-sum test)
        │     - Strict thresholding: |avg_log2FC| > 0.5, Bonferroni-adjusted p < 0.05
        │     - Custom GO database: eggNOG annotation + clusterProfiler::go2term
        │     - Background universe: all eggNOG-annotated genes with GO terms (BH FDR, Table S1/S4/S5)
        │
        ├─► 06_SCENIC_workflow.sh & 06_SCENIC_analysis.R
        │     - Cross-species orthology mapping with FlyBase D. melanogaster (DIAMOND, E < 1e-5; Table S2/S3)
        │     - pySCENIC (v0.12.1): GRNBoost2 -> RcisTarget -> AUCell regulon scoring
        │     - Integration of regulon AUCell scores into Seurat and stage-dynamic TF profiling
        │
        └─► 07_CellChat_communication.R
              - Cell-cell communication inference with CellChat (v1.6.1)
              - Custom crustacean-fly ligand-receptor database built from FlyPhoneDB2
              - Stage-resolved signaling networks (pre-, early-, and rapid-developmental stages, p < 0.05)
```

---

## Repository Structure

```text
.
├── .gitignore
├── README.md
├── scripts/
│   ├── 01_QC_and_harmony_integration.R
│   ├── 02_clustering_and_PC_sensitivity.R
│   ├── 03_EGC_subclustering_and_cellcycle.R
│   ├── 04_monocle3_trajectory.R
│   ├── 05_differential_expression_and_GO_KEGG.R
│   ├── 06_SCENIC_workflow.sh
│   ├── 06_SCENIC_analysis.R
│   └── 07_CellChat_communication.R
└── data/
    └── README.md (Instructions on downloading raw/processed data from public archives)
```

---

## System Requirements & Software Versions

The analysis scripts were developed and tested under Linux (Ubuntu 22.04 LTS / CentOS 7) and R (>= 4.3.0).

### Core Dependencies:
- **R packages**:
  - `Seurat` (v5.3.0)
  - `harmony` (v1.2.3)
  - `DoubletFinder` (v2.0.6)
  - `monocle3` (v1.4.27)
  - `clusterProfiler` (>= 4.10.0)
  - `CellChat` (v1.6.1)
  - `mclust` (>= 6.0.0)
  - `SCopeLoomR` (>= 0.13.0)
  - `AUCell` (>= 1.24.0)
  - `tidyverse` (`dplyr`, `ggplot2`, `readr`, `tibble`, `tidyr`)
  - `openxlsx` (>= 4.2.5)

- **Command-line tools & Python environment**:
  - `Cell Ranger` (v9.0.1)
  - `DIAMOND` (v2.1.8)
  - `pySCENIC` (v0.12.1)
  - Python (>= 3.9) with `loompy`, `scanpy`, `arboreto`

---

## Quick Start & Usage

### 1. Environment Setup
To reproduce the analysis environment using Conda:
```bash
# Create and activate R & Python environments
conda create -n es_snRNA python=3.9 r-base=4.3.0 -y
conda activate es_snRNA

# Install command-line tools
conda install -c bioconda diamond=2.1.8 -y
pip install pyscenic==0.12.1

# Install R packages within R
R -e "install.packages(c('Seurat', 'harmony', 'tidyverse', 'openxlsx', 'mclust', 'BiocManager'))"
R -e "BiocManager::install(c('clusterProfiler', 'GO.db', 'AUCell', 'SCopeLoomR', 'ComplexHeatmap'))"
R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', ref = 'v2.0.6')"
R -e "remotes::install_github('sqjin/CellChat', ref = 'v1.6.1')"
R -e "remotes::install_github('cole-trapnell-lab/monocle3', ref = 'v1.4.27')"
```

### 2. Execution Order
1. **Quality Control & Harmony Batch Integration**:
   ```bash
   Rscript scripts/01_QC_and_harmony_integration.R
   ```
2. **Major Cell Clustering & PC Sensitivity Evaluation (Figure S1)**:
   ```bash
   Rscript scripts/02_clustering_and_PC_sensitivity.R
   ```
3. **EGC Subclustering & Cell-Cycle Phase Assessment**:
   ```bash
   Rscript scripts/03_EGC_subclustering_and_cellcycle.R
   ```
4. **Pseudotime Trajectory Inference (Monocle3)**:
   ```bash
   Rscript scripts/04_monocle3_trajectory.R
   ```
5. **Differential Gene Expression & eggNOG-based Custom GO/KEGG Enrichment**:
   ```bash
   Rscript scripts/05_differential_expression_and_GO_KEGG.R
   ```
6. **Cross-species pySCENIC Regulatory Network Workflow**:
   ```bash
   # Run DIAMOND homology search and pySCENIC inference in bash
   bash scripts/06_SCENIC_workflow.sh
   # Downstream regulon activity integration and visualization in R
   Rscript scripts/06_SCENIC_analysis.R
   ```
7. **Cell-Cell Communication via FlyPhoneDB2 (CellChat)**:
   ```bash
   Rscript scripts/07_CellChat_communication.R
   ```

---

## Data Availability

- Raw sequencing reads (FASTQ) and processed expression matrices have been deposited in the NCBI GEO / SRA / GSA database under accession number `[GSEXXXXXX] / [PRJNAXXXXXX]`.
- Homology mapping tables (Table S2, Table S3) and full enrichment records (Table S1/S4/S5) are available as supplementary files accompanying the manuscript.

---

## Citation

If you use the scripts or resources in this repository, please cite:
> Yang et al., Single-nucleus transcriptomic profiling reveals cellular dynamics and regulatory networks during accessory sex gland development in the Chinese mitten crab (*Eriocheir sinensis*). (2025/2026, *Under Review*).
