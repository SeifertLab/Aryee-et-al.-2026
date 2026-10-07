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
| `orthology/Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv` | `512758AF98B6C20295B47AD59E04858918748874CBBB2DA2DC34D8CA4D9A3BF7` | Collaborator-supplied OrthoFinder pair table used by the paper |
| `metadata/sample_sheet_dual.csv` | `C53E593DC6D2D05E22C18D4BBC2AA131A8341E8E3E7AC43F1AEDFCAAD1CC9F5C` | Eighteen-sample species and age metadata |
| `figure_inputs/young_mito2_redox_panel_groups.csv` | `B3F5D7A0C094E21E71FBBFD1FA6A4CF6C32044CCF811ED28CBF99938BC09EA1F` | MitoCarta-derived pathway membership and gene order |

Hashes for repository text files are calculated with LF line endings, as specified in `.gitattributes`.

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
