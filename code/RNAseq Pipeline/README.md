# Raw-read processing commands

These scripts document the bulk-seq pipeline for rocessing raw fastqc files into the gene-count matrices used in the manuscript analysis:

1. `1_*_hisat2.sh` — align paired-end FASTQ files to the species reference with HISAT2.
2. `2_*_samtools.sh` — sort the SAM alignments into BAM files.
3. `3_*_stringtie_manual.sh` — generate a reference-guided GTF for each sample with StringTie and build the species merge list.
4. `3.1_*_stringtie_merge.sh` — merge the per-sample transcript annotations for the GffCompare QC branch.
5. `4_*_gffcompare.sh` — compare the merged transcript annotation with the species reference GTF. This is transcript-annotation QC and does not modify the deposited count matrices.
6. `5_*_transcript_abundance.sh` — quantify each sample with StringTie using the species genomic GTF for count extraction.
7. `6_download_prepde.sh` and `7_run_prepde.sh` — obtain StringTie's prepDE utility and create separate *Acomys* and *Mus* gene-count matrices.


The scripts retain the original Linux paths and require path edits before rerunning. Steps 3, 3.1, and 4 form a transcript-annotation QC branch, but steps 3.1 and 4 are not required for progression to step 5; step 5 supplies the per-sample GTF files used by `prepDE.py3`. 

The reproducible R workflow begins with the count matrices deposited in GEO GSE339424. 
