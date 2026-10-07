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

The shell scripts under `code/RNAseq Pipeline/` preserve the original paths and parameters used between alignment and count-matrix construction. FastQC and MultiQC commands are not included, and the scripts are not containerized. 

The supported reproducible workflow begins with the processed count matrices deposited in GEO GSE339424; all inputs needed for the manuscript R analysis after that point are identified in this repository.
