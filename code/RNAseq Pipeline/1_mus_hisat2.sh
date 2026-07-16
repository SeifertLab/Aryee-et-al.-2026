#!/bin/bash

# Define the path to the new directory containing the .fq files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/MUS"
OUTPUT_DIR="$DATA_DIR/output"

# Check if the new directory exists. If not, create it.
[ -d "$OUTPUT_DIR" ] || mkdir "$OUTPUT_DIR"

SCRIPT="Hisat2_align"
SCRIPT_NAME="Hisat2 align (-x) to paired end reads, _1.gz and _2.gz"
ERROR_LOG="$DATA_DIR/error_log_$SCRIPT.txt"

echo "Running $SCRIPT_NAME script on $DIR"
# Define a function to execute when there's an error, before the script exits.
trap 'echo "Error in $SCRIPT_NAME script at line $LINENO" >> $ERROR_LOG' ERR

# Define the path to the genome index
HISAT2_INDEX="/home/brennan/GRCm38/genome_index"

hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/MEF1_1.fq" -2 "$DATA_DIR/MEF1_2.fq" -S "$OUTPUT_DIR/aligned_MEF1.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/MEF2_1.fq" -2 "$DATA_DIR/MEF2_2.fq" -S "$OUTPUT_DIR/aligned_MEF2.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/MEF3_1.fq" -2 "$DATA_DIR/MEF3_2.fq" -S "$OUTPUT_DIR/aligned_MEF3.sam"

hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/MUSY1_1.fq" -2 "$DATA_DIR/MUSY1_2.fq" -S "$OUTPUT_DIR/aligned_MUSY1.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/MUSY2_1.fq" -2 "$DATA_DIR/MUSY2_2.fq" -S "$OUTPUT_DIR/aligned_MUSY2.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/MUSY3_1.fq" -2 "$DATA_DIR/MUSY3_2.fq" -S "$OUTPUT_DIR/aligned_MUSY3.sam"

hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/Mus_old1_1.fq" -2 "$DATA_DIR/Mus_old1_2.fq" -S "$OUTPUT_DIR/aligned_Mus_old1.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/Mus_Old2_1.fq" -2 "$DATA_DIR/Mus_Old2_2.fq" -S "$OUTPUT_DIR/aligned_Mus_Old2.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/Mus_Old3_1.fq" -2 "$DATA_DIR/Mus_Old3_2.fq" -S "$OUTPUT_DIR/aligned_Mus_Old3.sam"
