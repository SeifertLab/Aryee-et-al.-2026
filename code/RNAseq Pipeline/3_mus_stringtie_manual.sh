#!/bin/bash

# Define the path to the directory containing the .bam files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/MUS/output"
OUTPUT_DIR="/home/brennan/RNAseqAPRIL/clean_april/MUS/output"

# Define the path to the genome gtf
GTF="/home/brennan/GRCm38/Mus_musculus.GRCm38.102.gtf"


echo "Starting stringtie with whole directory, running manually per sample"

# Start a fresh merge list for step 3.1 (stringtie --merge)
> "$OUTPUT_DIR/lst_for_merge.txt"

stringtie -p 40 "$DATA_DIR/aligned_MEF1.bam" -G "$GTF" -e -l MEF1 -o "$OUTPUT_DIR/output_MEF1.gtf"
echo "$OUTPUT_DIR/output_MEF1.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "MEF1 done"
stringtie -p 40 "$DATA_DIR/aligned_MEF2.bam" -G "$GTF" -e -l MEF2 -o "$OUTPUT_DIR/output_MEF2.gtf"
echo "$OUTPUT_DIR/output_MEF2.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "MEF2 done"
stringtie -p 40 "$DATA_DIR/aligned_MEF3.bam" -G "$GTF" -e -l MEF3 -o "$OUTPUT_DIR/output_MEF3.gtf"
echo "$OUTPUT_DIR/output_MEF3.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "MEF3 done"

stringtie -p 40 "$DATA_DIR/aligned_MUSY1.bam" -G "$GTF" -e -l MUSY1 -o "$OUTPUT_DIR/output_MUSY1.gtf"
echo "$OUTPUT_DIR/output_MUSY1.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "MUSY1 done"
stringtie -p 40 "$DATA_DIR/aligned_MUSY2.bam" -G "$GTF" -e -l MUSY2 -o "$OUTPUT_DIR/output_MUSY2.gtf"
echo "$OUTPUT_DIR/output_MUSY2.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "MUSY2 done"
stringtie -p 40 "$DATA_DIR/aligned_MUSY3.bam" -G "$GTF" -e -l MUSY3 -o "$OUTPUT_DIR/output_MUSY3.gtf"
echo "$OUTPUT_DIR/output_MUSY3.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "MUSY3 done"


stringtie -p 40 "$DATA_DIR/aligned_Mus_old1.bam" -G "$GTF" -e -l Mus_old1 -o "$OUTPUT_DIR/output_Mus_old1.gtf"
echo "$OUTPUT_DIR/output_Mus_old1.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "Mus_old1 done"
stringtie -p 40 "$DATA_DIR/aligned_Mus_Old2.bam" -G "$GTF" -e -l Mus_Old2 -o "$OUTPUT_DIR/output_Mus_Old2.gtf"
echo "$OUTPUT_DIR/output_Mus_Old2.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "Mus_Old2 done"
stringtie -p 40 "$DATA_DIR/aligned_Mus_Old3.bam" -G "$GTF" -e -l Mus_Old3 -o "$OUTPUT_DIR/output_Mus_Old3.gtf"
echo "$OUTPUT_DIR/output_Mus_Old3.gtf" >> "$OUTPUT_DIR/lst_for_merge.txt"
echo "Mus_Old3 done"
