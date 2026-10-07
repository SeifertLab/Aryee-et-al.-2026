# Script versions

Release candidate: `v1.0.0-rc2`
Register date: 2026-10-07
Starting public commit: `6452a87646b57f56286c22f794d8d2275bad8cd7`

SHA-256 hashes are calculated with LF line endings, as specified in `.gitattributes`.

## Manuscript R workflow

| Script | SHA-256 | Role |
|---|---|---|
| `code/R/run_young_only_dual_species_deseq2.R` | `FBC7A83D35656460BED72472AF7095920D01560E8274743D16A3A089DCDFFBC7` | Count alignment, strict reciprocal 1:1 ortholog filtering, count filtering, young-only DESeq2 |
| `code/R/render_rotated_sig_blocks_inferno.R` | `C269546C2C519DAED668A9C1B4A68AF68E330B841288064E7F64C90ADA7598C2` | Main five-block significant-only heatmap |
| `code/R/render_standalone_fission_fusion_reviewer_versions.R` | `1A358737CC60B1558EC9E03F3192F4DE4F2B4F7A84464D8ECA9780AFE76DE680` | Supplemental fission/fusion heatmaps |

## DE analysis inputs

| Input | SHA-256 |
|---|---|
| `data/orthology/Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv` | `512758AF98B6C20295B47AD59E04858918748874CBBB2DA2DC34D8CA4D9A3BF7` |
| `data/metadata/sample_sheet_dual.csv` | `C53E593DC6D2D05E22C18D4BBC2AA131A8341E8E3E7AC43F1AEDFCAAD1CC9F5C` |
| `data/figure_inputs/young_mito2_redox_panel_groups.csv` | `B3F5D7A0C094E21E71FBBFD1FA6A4CF6C32044CCF811ED28CBF99938BC09EA1F` |

The GEO count-matrix hashes are recorded in `data/README.md`.

## Bulk-seq pipeline

These files record the commands used between FASTQ alignment and count-matrix construction. They retain the original machine paths and are not a portable workflow. FastQC and MultiQC commands are not included. The deposited GEO count matrices are the supported input for the R workflow.

| Script | SHA-256 | Role |
|---|---|---|
| `code/RNAseq Pipeline/1_aco_hisat2.sh` | `C9A8DD46C38A9F0F4D61FBED3017A16EC3C37644FDEA070067AA7FBCFB28DB85` | Acomys paired-end alignment |
| `code/RNAseq Pipeline/1_mus_hisat2.sh` | `53715679FDD3C26912E79246E82D6D57724A7FAD56A6845CA2FCA00937351A68` | Mus paired-end alignment |
| `code/RNAseq Pipeline/2_aco_samtools.sh` | `B16A8C7AF3A1842562EC8652D969ED5D23353241479B5D75D940E8F74DF9EF16` | Acomys SAM-to-sorted-BAM conversion |
| `code/RNAseq Pipeline/2_mus_samtools.sh` | `A34B4DEC0D8394116F4CB6AEFBB3C7664A2AEB99A12F2B4E276A1E41E53FEAEC` | Mus SAM-to-sorted-BAM conversion |
| `code/RNAseq Pipeline/3_aco_stringtie_manual.sh` | `78E32EF5F8DE225413977D846D963127278FD0EBD204EFCB759DDD439AC90DB8` | Acomys per-sample reference-guided GTF generation and merge-list construction |
| `code/RNAseq Pipeline/3_mus_stringtie_manual.sh` | `C427F68B679060A9E907E54C99C3CEA874A179AE158F90F709F002DCB87FE567` | Mus per-sample reference-guided GTF generation and merge-list construction |
| `code/RNAseq Pipeline/3.1_aco_stringtie_merge.sh` | `8A290714B10F349889632F2C3481416F78CCB5732ECE5BA30E9876ED576D76C0` | Acomys transcript-annotation merge for GffCompare QC |
| `code/RNAseq Pipeline/3.1_mus_stringtie_merge.sh` | `EFC79BD7D4EA7B22798AB3BD1634D43ED773FE1619E85EE97B5C8D926A31F2B4` | Mus transcript-annotation merge for GffCompare QC |
| `code/RNAseq Pipeline/4_aco_gffcompare.sh` | `2F10A929A0EAA8DAE6DB9E422A6C51A2624A9FE050F7A4924CD68CB5F5CEB61E` | Acomys merged-annotation comparison against the reference GTF |
| `code/RNAseq Pipeline/4_mus_gffcompare.sh` | `3A166C8F3212566D64C4B043297FF9303D25E1887FEE4B9B96DD83907EA04745` | Mus merged-annotation comparison against the reference GTF |
| `code/RNAseq Pipeline/5_aco_transcript_abundance.sh` | `CA4AEEFD032E80F106372C6E19DC34762421FDBC2274E7F885EEC94C21449AE8` | Acomys reference-guided StringTie quantification |
| `code/RNAseq Pipeline/5_mus_transcript_abundance.sh` | `2E10FD3179433F9EC4BAB9D086786F73A0E12B550F046522F6F239BE2510B9B0` | Mus reference-guided StringTie quantification |
| `code/RNAseq Pipeline/6_download_prepde.sh` | `46A2719D880E06AFE91E11CE5396B6DFDA95ACBBAFCEF35DD758B1B6C82A20A2` | prepDE download command |
| `code/RNAseq Pipeline/7_run_prepde.sh` | `E918929F9D402470502669C0F9FCB7C6589701A94B4D00C35A43E932CF7F237F` | Builds the separate Acomys and Mus gene-count matrices |
