#!/bin/bash

# Define the path to the new directory containing the .fq files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO"
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
HISAT2_INDEX="/home/brennan/ACO_GCA_907164435.1/genome_index"

# i.e. every *_1.fq / *_2.fq pair in DATA_DIR: the 9 Acomys libraries.
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/AEF1_1.fq" -2 "$DATA_DIR/AEF1_2.fq" -S "$OUTPUT_DIR/aligned_AEF1.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/AEF2_1.fq" -2 "$DATA_DIR/AEF2_2.fq" -S "$OUTPUT_DIR/aligned_AEF2.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/AEF3_1.fq" -2 "$DATA_DIR/AEF3_2.fq" -S "$OUTPUT_DIR/aligned_AEF3.sam"

hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/Aco_5515_1.fq" -2 "$DATA_DIR/Aco_5515_2.fq" -S "$OUTPUT_DIR/aligned_Aco_5515.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/Aco_5516_1.fq" -2 "$DATA_DIR/Aco_5516_2.fq" -S "$OUTPUT_DIR/aligned_Aco_5516.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/Aco_5517_1.fq" -2 "$DATA_DIR/Aco_5517_2.fq" -S "$OUTPUT_DIR/aligned_Aco_5517.sam"

hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/SMA1_1.fq" -2 "$DATA_DIR/SMA1_2.fq" -S "$OUTPUT_DIR/aligned_SMA1.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/SMA5_1.fq" -2 "$DATA_DIR/SMA5_2.fq" -S "$OUTPUT_DIR/aligned_SMA5.sam"
hisat2 -p 36 --dta -x "$HISAT2_INDEX" -1 "$DATA_DIR/SMA6_1.fq" -2 "$DATA_DIR/SMA6_2.fq" -S "$OUTPUT_DIR/aligned_SMA6.sam"
