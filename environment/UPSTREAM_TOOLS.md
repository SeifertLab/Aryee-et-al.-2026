# Upstream RNA-seq processing

The manuscript used the following raw-read processing:

| Stage | Tool/version |
|---|---|
| Read quality control | FastQC 0.11.9 |
| QC aggregation | MultiQC 1.16 |
| Alignment | HISAT2 2.2.1 |
| SAM/BAM processing | SAMtools 1.6 |
| Reference-guided transcript quantification | StringTie 2.1.7 |
| Gene-count matrix construction | StringTie `prepDE.py3` |

References were *Acomys dimidiatus* Ensembl Rapid Release assembly GCA_907164435.1 and *Mus musculus* GRCm38 with Ensembl release 102 gene models. The analysis used annotation-guided per-sample StringTie estimates.

The bulk RNAseq pipeline shell commands under `code/RNAseq Pipeline/` preserve the original machine paths and parameters used to generate the count matrices. They are not containerized. Consequently, the public entry point is the pair of processed count matrices deposited in GEO GSE339424. The manuscript R workflow is fully defined from those matrices.