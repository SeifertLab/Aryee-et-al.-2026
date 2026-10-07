# Aryee et al. RNA-seq analysis and figure code

This repository contains the RNA-seq differential-expression analysis and the two heatmap renderers used for the published article:

> Aryee E, Aloysius A, Saxena S, Donahue R, Riddell B, Vekaria H, Sullivan PG, Liu H, Kachroo P, Seifert AW. *Spiny mouse fibroblasts exhibit a baseline preference for glycolysis and are resilient to oxidative stress across lifespan.* Journal of Biological Chemistry (2026). [Published article](https://www.sciencedirect.com/science/article/pii/S0021925826024774), [DOI](https://doi.org/10.1016/j.jbc.2026.113605)

The workflow starts from the species-specific gene-count matrices deposited in NCBI GEO [GSE339424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE339424). It constructs the cross-species gene universe from the OrthoFinder alignment used in the paper, runs the young *Mus musculus* versus young *Acomys dimidiatus* DESeq2 comparison, and renders the manuscript heatmaps.

## Analysis and figure scripts

Run these three scripts in order:

1. `code/R/run_young_only_dual_species_deseq2.R` — constructs the strict reciprocal 1:1 ortholog count matrix and performs the six-sample young-only DESeq2 analysis.
2. `code/R/render_rotated_sig_blocks_inferno.R` — renders the main significant-only metabolic/redox heatmap.
3. `code/R/render_standalone_fission_fusion_reviewer_versions.R` — renders the supplemental mitochondrial fission/fusion heatmaps.

## Inputs

The orthology table, sample metadata, and figure-panel definitions used by the analysis are included in this repository:

- `data/orthology/Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv` — pairwise OrthoFinder output used to define strict reciprocal 1:1 orthologs;
- `data/metadata/sample_sheet_dual.csv` — metadata: sample identities, species, and age groups;
- `data/figure_inputs/young_mito2_redox_panel_groups.csv` — curated, MitoCarta-derived pathway membership and gene order used in the heatmaps.

Download the two processed count matrices from GEO and place them in `data/processed_counts/` as described in [data/README.md](data/README.md).

## Requirements

The analysis used R 4.5.2 with:

| Package | Version |
|---|---:|
| DESeq2 | 1.50.2 |
| apeglm | 1.32.0 |
| dplyr | 1.2.0 |
| readr | 2.2.0 |
| stringr | 1.6.0 |
| tidyr | 1.3.2 |
| tibble | 3.3.1 |
| ComplexHeatmap | 2.26.1 |
| circlize | 0.4.18 |

Install the required packages with:

```r
install.packages(c("BiocManager", "dplyr", "readr", "stringr", "tidyr", "tibble", "circlize"))
BiocManager::install(c("DESeq2", "apeglm", "ComplexHeatmap"))
```

The R workflow was validated on Windows 11 with R 4.5.2. It has not been validated in other R versions or operating systems.

## Reproduce the analysis

```bash
git clone https://github.com/SeifertLab/Aryee-et-al.-2026.git
cd Aryee-et-al.-2026

# Add the two GEO count matrices under data/processed_counts/, then run:
Rscript code/R/run_young_only_dual_species_deseq2.R
Rscript code/R/render_rotated_sig_blocks_inferno.R
Rscript code/R/render_standalone_fission_fusion_reviewer_versions.R
```

Default outputs are written to:

- `results/young_only_deseq2/`
- `figures/main/`
- `figures/supplement/`

See [docs/RUNBOOK.md](docs/RUNBOOK.md) for the expected files and numerical checks.

## Repository structure

```text
code/
├── R/                    # DESeq2 analysis and manuscript figure renderers
└── RNAseq Pipeline/      # FASTQ-to-count pipeline commands
data/
├── figure_inputs/        # Curated figure membership and gene order
├── metadata/             # Sample identities, species, and age groups
├── orthology/            # OrthoFinder mapping and provenance
└── processed_counts/     # Place the two GEO count matrices here
docs/                     # Runbook, version register, and validation record
environment/              # R and upstream tool versions
```

## Statistical analysis

The OrthoFinder table is expanded to individual gene pairs after removing Ensembl version suffixes. A pair is retained only when the *Acomys* gene has one *Mus* partner and the *Mus* gene has one *Acomys* partner. Of 16,314 strict reciprocal pairs, 16,310 are represented in both prepDE count matrices. The paper's filter retains genes with at least 10 counts in at least 3 of all 18 samples, yielding 12,944 tested ortholog pairs.

The primary contrast contains three *Acomys* samples (SMA1, SMA5, SMA6) and three *Mus* samples (MUSY1, MUSY2, MUSY3). DESeq2 uses the design `~ species`, with *Acomys* as the reference - positive log2 fold change therefore means higher expression in *Mus*. Wald-test p-values are adjusted by the Benjamini-Hochberg method and effect sizes are shrunk with apeglm.

Heatmaps show row-wise z scores of `log2(DESeq2 normalized count + 1)` across the six samples. The main heatmap includes genes from the curated panel with BH FDR < 0.05.

## Validation

The validated release produces 12,944 tested ortholog pairs, of which 9,017 had BH FDR < 0.05. Details are in [docs/VALIDATION.md](docs/VALIDATION.md).

## Upstream sequencing workflow

The deposited count matrices were generated from paired-end reads using FastQC (0.11.9), MultiQC (1.16), HISAT2 (2.2.1), SAMtools (1.6), StringTie (2.1.7), and StringTie's `prepDE.py3`, with *Acomys* assembly GCA_907164435.1 and *Mus* GRCm38/Ensembl release 102 annotations. The original processing commands are retained under `code/RNAseq Pipeline/`. Those scripts contain machine-specific paths and are not portable without modifying these for your device.

## Data availability

- Sequencing data and processed count matrices: NCBI GEO [GSE339424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE339424).
- Published article: [ScienceDirect/JBC](https://www.sciencedirect.com/science/article/pii/S0021925826024774).

Source code is maintained at [SeifertLab/Aryee-et-al.-2026](https://github.com/SeifertLab/Aryee-et-al.-2026). A versioned archival DOI will be added here after the GitHub release is deposited in Zenodo.

## Citation

PLACEHOLDER FOR ZENODO DOI

## Acknowledgements

This research was funded in part by Aligning Science Across Parkinson’s (ASAP-020495) through The Michael J. Fox Foundation for Parkinson’s Research (MJFF).

Brennan Riddell prepared the public analysis workflow and repository. We thank Huayun Chen for generating and providing the *Acomys dimidiatus–Mus musculus* OrthoFinder table used to construct the cross-species gene universe.

## License

This code is released under the MIT license.
