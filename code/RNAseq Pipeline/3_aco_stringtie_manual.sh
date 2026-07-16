#!/bin/bash

# Define the path to the DATA_DIRectory containing the .bam files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO/output"
OUTPUT_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO/output"

# Define the path to the genome gtf
GTF="/home/brennan/ACO_GCA_907164435.1/Acomys_dimidiatus-GCA_907164435.1-2022_07-genes.gtf"

echo "Starting stringtie with whole directory, running manually per sample"

# Start a fresh merge list for step 3.1 (stringtie --merge)
> "$OUTPUT_DIR/lst_for_merge.txt"

stringtie -o "$OUTPUT_DIR/output_AEF1.gtf" -p 40 -G "$GTF" -e -l AEF1 "$DATA_DIR/aligned_AEF1.bam"
echo "$OUTPUT_DIR/output_AEF1.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "AEF1 done"
stringtie -o "$OUTPUT_DIR/output_AEF2.gtf" -p 40 -G "$GTF" -e -l AEF2 "$DATA_DIR/aligned_AEF2.bam"
echo "$OUTPUT_DIR/output_AEF2.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "AEF2 done"
stringtie -o "$OUTPUT_DIR/output_AEF3.gtf" -p 40 -G "$GTF" -e -l AEF3 "$DATA_DIR/aligned_AEF3.bam"
echo "$OUTPUT_DIR/output_AEF3.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "AEF3 done"

stringtie -o "$OUTPUT_DIR/output_SMA1.gtf" -p 40 -G "$GTF" -e -l SMA1 "$DATA_DIR/aligned_SMA1.bam"
echo "$OUTPUT_DIR/output_SMA1.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "SMA1 done"
stringtie -o "$OUTPUT_DIR/output_SMA5.gtf" -p 40 -G "$GTF" -e -l SMA5 "$DATA_DIR/aligned_SMA5.bam"
echo "$OUTPUT_DIR/output_SMA5.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "SMA5 done"
stringtie -o "$OUTPUT_DIR/output_SMA6.gtf" -p 40 -G "$GTF" -e -l SMA6 "$DATA_DIR/aligned_SMA6.bam"
echo "$OUTPUT_DIR/output_SMA6.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "SMA6 done"


stringtie -o "$OUTPUT_DIR/output_Aco_5515.gtf" -p 40 -G "$GTF" -e -l Aco_5515 "$DATA_DIR/aligned_Aco_5515.bam"
echo "$OUTPUT_DIR/output_Aco_5515.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "Aco_5515 done"
stringtie -o "$OUTPUT_DIR/output_Aco_5516.gtf" -p 40 -G "$GTF" -e -l Aco_5516 "$DATA_DIR/aligned_Aco_5516.bam"
echo "$OUTPUT_DIR/output_Aco_5516.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "Aco_5516 done"
stringtie -o "$OUTPUT_DIR/output_Aco_5517.gtf" -p 40 -G "$GTF" -e -l Aco_5517 "$DATA_DIR/aligned_Aco_5517.bam"
echo "$OUTPUT_DIR/output_Aco_5517.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "Aco_5517 done"
