#!/bin/bash

# Path to the directory with the generated files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO/output"
OUT_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO/output"
GTF="/home/brennan/ACO_GCA_907164435.1/Acomys_dimidiatus-GCA_907164435.1-2022_07-genes.gtf"

stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/AEF1_genome_abundance.gtf" "$DATA_DIR/aligned_AEF1.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/AEF2_genome_abundance.gtf" "$DATA_DIR/aligned_AEF2.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/AEF3_genome_abundance.gtf" "$DATA_DIR/aligned_AEF3.bam"

stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/SMA1_genome_abundance.gtf" "$DATA_DIR/aligned_SMA1.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/SMA5_genome_abundance.gtf" "$DATA_DIR/aligned_SMA5.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/SMA6_genome_abundance.gtf" "$DATA_DIR/aligned_SMA6.bam"

stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/Aco_5515_genome_abundance.gtf" "$DATA_DIR/aligned_Aco_5515.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/Aco_5516_genome_abundance.gtf" "$DATA_DIR/aligned_Aco_5516.bam"
stringtie -e -B -p 40 -G "$GTF" -o "$OUT_DIR/Aco_5517_genome_abundance.gtf" "$DATA_DIR/aligned_Aco_5517.bam"
