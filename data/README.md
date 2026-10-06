# Analysis inputs

## Download from GEO

Download the processed prepDE gene-count matrices from NCBI GEO accession [GSE339424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE339424) and place them under `data/processed_counts/` with these names:

| File | SHA-256 |
|---|---|
| `gene_count_matrix_GTF_acomys.csv` | `7F6F1DA89BE3F567683899C0655D01F7A6FCDF032A9276B8F9A1A1AA0694F2C9` |
| `gene_count_matrix_GTF_mus.csv` | `CC13CEE214C3FFEF2762333D2017685FB6E2EFC0529E36D713F3C630E3E8ABBA` |

The count matrices are intentionally not duplicated in Git.

## Versioned with the code

| File | SHA-256 | Role |
|---|---|---|
| `orthology/Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv` | `512758AF98B6C20295B47AD59E04858918748874CBBB2DA2DC34D8CA4D9A3BF7` | Exact collaborator-supplied OrthoFinder pair table used by the paper |
| `metadata/sample_sheet_dual.csv` | `0CC1E52E622BDD026E2517497B60666790A5C8D2B2FCC17D4914E9286D1E50EF` | Eighteen-sample species and age metadata |
| `figure_inputs/young_mito2_redox_panel_groups.csv` | `ABD49A085AA19B3B8C94795BD62788FAA09599E0F8D3FDD04D869CC06717815F` | Frozen first-author-curated, MitoCarta-derived figure panel |

## Expected layout

```text
data/
├── figure_inputs/
│   └── young_mito2_redox_panel_groups.csv
├── metadata/
│   └── sample_sheet_dual.csv
├── orthology/
│   ├── Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv
│   └── ORTHOLOGY_PROVENANCE.md
└── processed_counts/
    ├── README.md
    ├── gene_count_matrix_GTF_acomys.csv
    └── gene_count_matrix_GTF_mus.csv
```
