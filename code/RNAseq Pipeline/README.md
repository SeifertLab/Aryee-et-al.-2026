# Historical raw-read processing commands

These scripts document the direct processing path used to create the GEO gene-count matrices:

1. `1_*_hisat2.sh` — align paired-end FASTQ files to the species reference with HISAT2.
2. `2_*_samtools.sh` — sort the SAM alignments into BAM files.
3. `3_*_stringtie_manual.sh` - reference guided alignment of BAM files to genomic .gtf
4. `4_*_gffcompare.sh` - GFFcompare, QC step for alignment accuracy - does not change output
5. `5_*_transcript_abundance.sh` — quantify each sample with StringTie using the species genomic GTF.
6. `6_download_prepde.sh` and `7_run_prepde.sh` — obtain StringTie's prepDE utility and create separate *Acomys* and *Mus* gene-count matrices.


The scripts retain the original Linux paths and should be treated as provenance. Edit the path variables before rerunning. The supported public downstream workflow begins with the count matrices deposited in GEO GSE339424.
