# Define the path to the directory containing the .sam files
DATA_DIR="/home/brennan/RNAseqAPRIL/clean_april/ACO/output"

SCRIPT="Samtools"
SCRIPT_NAME="Samtools - .sam to .bam"
ERROR_LOG="$DATA_DIR/error_log_$SCRIPT.txt"
echo "Running $SCRIPT_NAME script on $DATA_DIR"
# Define a function to execute when there's an error, before the script exits.
trap 'echo "Error in $SCRIPT_NAME script at line $LINENO" >> $ERROR_LOG' ERR

# Iterate over .sam files in the directory
for SAM_FILE in "$DATA_DIR"/*.sam; do
    # Get the base name of the .sam file without path and extension
    BASENAME=$(basename "$SAM_FILE" .sam)

    # Run the samtools command
    samtools sort -@ 8 -o "$DATA_DIR/$BASENAME.bam" "$SAM_FILE"
done
