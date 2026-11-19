# Data Directory Structure and Description

This folder should contain or link to the input datasets required to reproduce the analysis. Due to file size limits, large raw files (.fastq.gz, .bam) and intermediate binary objects (.rds, .loom) are not hosted directly in this GitHub repository.

## Recommended Directory Structure

```text
data/
├── raw/
│   ├── sample_ASG1/filtered_feature_bc_matrix/
│   ├── sample_ASG2/filtered_feature_bc_matrix/
│   └── sample_ASG3/filtered_feature_bc_matrix/
├── reference/
│   ├── E_sinensis_genome.fasta
│   ├── E_sinensis_genes.gtf
│   ├── dmel-all-translation-r6.50.fasta (FlyBase Drosophila proteins)
│   └── esfemale.out.emapper.annotations.csv (eggNOG functional annotations)
└── metadata/
    ├── Table_S2_crab_fly_homology.csv
    └── Table_S3_TF_tier_classification.xlsx
```

## Public Data Accession
- Gene Expression Omnibus (GEO): `GSEXXXXXX`
- Genome Sequence Archive (GSA): `CRAXXXXXX`
- E. sinensis Reference Genome: GenBank `GCA_013436485.1`
