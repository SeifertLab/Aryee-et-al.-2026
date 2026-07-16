#!/bin/bash

# Define the path to the directory containing the .bam files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/MUS/output"

# Define the path to the genome gtf
GTF="/home/brennan/GRCm38/Mus_musculus.GRCm38.102.gtf"

# Run the stringtie merge command
stringtie --merge -p 40 -G "$GTF" -o "$DATA_DIR/stringtie_merged.gtf" "$DATA_DIR/lst_for_merge.txt"

