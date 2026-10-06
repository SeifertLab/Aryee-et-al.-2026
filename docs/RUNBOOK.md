# Reproduction runbook

## 1. Obtain the two GEO count matrices

Download the processed Acomys and Mus prepDE gene-count matrices from GEO GSE339424, rename them as specified in `data/README.md`, and place them in `data/processed_counts/`.

Verify their SHA-256 hashes before analysis:

In powershell:
```powershell
Get-FileHash data/processed_counts/gene_count_matrix_GTF_acomys.csv -Algorithm SHA256
Get-FileHash data/processed_counts/gene_count_matrix_GTF_mus.csv -Algorithm SHA256
```

## 2. Run the DE analysis

```bash
Rscript code/R/run_young_only_dual_species_deseq2.R
```

The script enforces these quality gates:

- 16,314 strict reciprocal 1:1 pairs in the frozen OrthoFinder table;
- 16,310 pairs represented in both species' count matrices;
- 12,944 pairs after `count >= 10 in at least 3 of 18 samples`;
- six young samples, three per species;
- one unambiguous `Mus versus Acomys` DESeq2 coefficient.

Expected files under `results/young_only_deseq2/`:

- `analysis_summary.csv`
- `ortholog_pair_presence_audit.csv`
- `tested_ortholog_manifest.csv`
- `aligned_filtered_raw_counts.tsv.gz`
- `young_only_dual_species_deseq2.csv`
- `young_only_normalized_counts.tsv.gz`
- `young_only_dds.rds`
- `young_only_results_raw.rds`
- `young_only_results_shrunk.rds`

## 3. Render the main figure heatmap

```bash
Rscript code/R/render_rotated_sig_blocks_inferno.R
```

This writes PDF, SVG, PNG, displayed-gene, and inclusion-audit files under `figures/main/`.

## 4. Render supplemental fission/fusion heatmaps

```bash
Rscript code/R/render_standalone_fission_fusion_reviewer_versions.R
```

This writes two figure variants and their audit tables under `figures/supplement/`:

- significant genes plus reviewer-requested Opa1, Dnm1l, Oma1, and Gdap1;
- every curated fission/fusion gene, with significance denoted.