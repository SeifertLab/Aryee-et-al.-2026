#!/bin/bash

SCRIPT_NAME="Download_prepDE.py3"
DIR="/home/brennan/RNAseqPLEASE/"  # current directory
ERROR_LOG="$DIR/error_log_$SCRIPT_NAME.txt"

echo "Running $SCRIPT_NAME script on $DIR"

# Define a function to execute when there's an error, before the script exits.
trap 'echo "Error in $SCRIPT_NAME script at line $LINENO" >> $ERROR_LOG' ERR

# The main logic of your script
curl -O https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 || exit 1
chmod +x prepDE.py3

echo "$SCRIPT_NAME completed successfully"
