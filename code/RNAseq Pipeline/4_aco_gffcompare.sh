#!/bin/bash

# Path to the directory with the merged file
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO/output"

# Path to the original genome annotation
GTF="/home/brennan/ACO_GCA_907164435.1/Acomys_dimidiatus-GCA_907164435.1-2022_07-genes.gtf"

SCRIPT="7.GFFCompare"
SCRIPT_NAME="7.GFFCompare"
ERROR_LOG="$DATA_DIR/error_log_$SCRIPT.txt"
echo "Running $SCRIPT_NAME script on $DATA_DIR"

# Define a function to execute when there's an error, before the script exits.
trap 'echo "Error in $SCRIPT_NAME script at line $LINENO" >> $ERROR_LOG' ERR

# Run gffcompare
#gffcompare -r "$GTF" -G -o "$DATA_DIR/merged" "$DATA_DIR/stringtie_merged.gtf"
gffcompare -R -r "$GTF" -o "$DATA_DIR/merged" "$DATA_DIR/stringtie_merged.gtf"
