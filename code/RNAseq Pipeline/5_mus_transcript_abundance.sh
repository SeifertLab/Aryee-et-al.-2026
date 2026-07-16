#!/bin/bash

# Path to the directory with the generated files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/MUS/output"
OUT_DIR="/home/brennan/RNAseqAPRIL/clean_april/MUS/output"
GTF="/home/brennan/GRCm38/Mus_musculus.GRCm38.102.gtf"

stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/MUSY1_genome_abundance.gtf" "$DATA_DIR/aligned_MUSY1.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/MUSY2_genome_abundance.gtf" "$DATA_DIR/aligned_MUSY2.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/MUSY3_genome_abundance.gtf" "$DATA_DIR/aligned_MUSY3.bam"

stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/MEF1_genome_abundance.gtf" "$DATA_DIR/aligned_MEF1.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/MEF2_genome_abundance.gtf" "$DATA_DIR/aligned_MEF2.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/MEF3_genome_abundance.gtf" "$DATA_DIR/aligned_MEF3.bam"

stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/Mus_old1_genome_abundance.gtf" "$DATA_DIR/aligned_Mus_old1.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/Mus_Old2_genome_abundance.gtf" "$DATA_DIR/aligned_Mus_Old2.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/Mus_Old3_genome_abundance.gtf" "$DATA_DIR/aligned_Mus_Old3.bam"
