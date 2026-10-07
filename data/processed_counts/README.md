# Processed gene-count matrices

Download the processed prepDE gene-count matrices from NCBI GEO accession GSE339424 and place them in this directory with these exact names:

| File | SHA-256 |
|---|---|
| `gene_count_matrix_GTF_acomys.csv` | `7F6F1DA89BE3F567683899C0655D01F7A6FCDF032A9276B8F9A1A1AA0694F2C9` |
| `gene_count_matrix_GTF_mus.csv` | `CC13CEE214C3FFEF2762333D2017685FB6E2EFC0529E36D713F3C630E3E8ABBA` |

These matrices were produced with StringTie quantification against the species' genomic GTF followed by `prepDE.py3`. They contain raw integer gene counts and are the entry point for the downstream R analysis.

