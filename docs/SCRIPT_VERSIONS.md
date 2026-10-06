# Script versions

Release candidate: `v1.0.0-rc2`  
Register date: 2026-10-06  
Starting public commit: `6452a87646b57f56286c22f794d8d2275bad8cd7`

## Manuscript R workflow

| Script | SHA-256 | Role |
|---|---|---|
| `code/R/run_young_only_dual_species_deseq2.R` | `99783D4725C6E10F13D7419797CC313C323B01BAC4AB9F582B2CC55DAC26A625` | Count alignment, strict reciprocal 1:1 ortholog filtering, count filtering, young-only DESeq2 |
| `code/R/render_rotated_sig_blocks_inferno.R` | `C269546C2C519DAED668A9C1B4A68AF68E330B841288064E7F64C90ADA7598C2` | Main five-block significant-only heatmap |
| `code/R/render_standalone_fission_fusion_reviewer_versions.R` | `ED39CCB438C66F568E7A17499EF8277B423D116C8961E16BBA2676C913834E35` | Supplemental fission/fusion heatmaps |

## DE analysis inputs

| Input | SHA-256 |
|---|---|
| `data/orthology/Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv` | `512758AF98B6C20295B47AD59E04858918748874CBBB2DA2DC34D8CA4D9A3BF7` |
| `data/metadata/sample_sheet_dual.csv` | `0CC1E52E622BDD026E2517497B60666790A5C8D2B2FCC17D4914E9286D1E50EF` |
| `data/figure_inputs/young_mito2_redox_panel_groups.csv` | `ABD49A085AA19B3B8C94795BD62788FAA09599E0F8D3FDD04D869CC06717815F` |

The GEO count-matrix hashes are recorded in `data/README.md`.

## Bulk-seq pipeline

These commands constitute the pipeline from raw FASTQ files to the deposited processed count matrices. They retain the original machine paths and are mirrored in the methodology in our NCBI GEO accession.

| Script | SHA-256 | Role |
|---|---|---|
| `code/RNAseq Pipeline/1_aco_hisat2.sh` | `C9A8DD46C38A9F0F4D61FBED3017A16EC3C37644FDEA070067AA7FBCFB28DB85` | Acomys paired-end alignment |
| `code/RNAseq Pipeline/1_mus_hisat2.sh` | `53715679FDD3C26912E79246E82D6D57724A7FAD56A6845CA2FCA00937351A68` | Mus paired-end alignment |
| `code/RNAseq Pipeline/2_aco_samtools.sh` | `B16A8C7AF3A1842562EC8652D969ED5D23353241479B5D75D940E8F74DF9EF16` | Acomys SAM-to-sorted-BAM conversion |
| `code/RNAseq Pipeline/2_mus_samtools.sh` | `A34B4DEC0D8394116F4CB6AEFBB3C7664A2AEB99A12F2B4E276A1E41E53FEAEC` | Mus SAM-to-sorted-BAM conversion |
| `code/RNAseq Pipeline/5_aco_transcript_abundance.sh` | `CA4AEEFD032E80F106372C6E19DC34762421FDBC2274E7F885EEC94C21449AE8` | Acomys reference-guided StringTie quantification |
| `code/RNAseq Pipeline/5_mus_transcript_abundance.sh` | `2E10FD3179433F9EC4BAB9D086786F73A0E12B550F046522F6F239BE2510B9B0` | Mus reference-guided StringTie quantification |
| `code/RNAseq Pipeline/6_download_prepde.sh` | `46A2719D880E06AFE91E11CE5396B6DFDA95ACBBAFCEF35DD758B1B6C82A20A2` | prepDE download command |
| `code/RNAseq Pipeline/7_run_prepde.sh` | `E918929F9D402470502669C0F9FCB7C6589701A94B4D00C35A43E932CF7F237F` | Builds the separate Acomys and Mus gene-count matrices |

