#!/usr/bin/env bash
set -euo pipefail

# Build species-specific count matrices from the StringTie GTF files in step 5.
PREPDE=${PREPDE:-./prepDE.py3}
ACO_GTF_DIR=${ACO_GTF_DIR:-/home/brennan/RNAseqAPRIL/clean_april/ACO/output}
MUS_GTF_DIR=${MUS_GTF_DIR:-/home/brennan/RNAseqAPRIL/clean_april/MUS/output}
COUNT_OUT_DIR=${COUNT_OUT_DIR:-.}

if [[ ! -f "$PREPDE" ]]; then
  echo "prepDE.py3 not found: $PREPDE" >&2
  exit 1
fi

aco_list=$(mktemp)
mus_list=$(mktemp)
trap 'rm -f "$aco_list" "$mus_list"' EXIT

for sample in AEF1 AEF2 AEF3 SMA1 SMA5 SMA6 Aco_5515 Aco_5516 Aco_5517; do
  printf '%s,%s/%s_genome_abundance.gtf\n' "$sample" "$ACO_GTF_DIR" "$sample" >> "$aco_list"
done

for sample in MEF1 MEF2 MEF3 MUSY1 MUSY2 MUSY3 Mus_old1 Mus_Old2 Mus_Old3; do
  printf '%s,%s/%s_genome_abundance.gtf\n' "$sample" "$MUS_GTF_DIR" "$sample" >> "$mus_list"
done

mkdir -p "$COUNT_OUT_DIR"

python3 "$PREPDE" \
  -i "$aco_list" \
  -g "$COUNT_OUT_DIR/gene_count_matrix_GTF_acomys.csv" \
  -t "$COUNT_OUT_DIR/transcript_count_matrix_GTF_acomys.csv"

python3 "$PREPDE" \
  -i "$mus_list" \
  -g "$COUNT_OUT_DIR/gene_count_matrix_GTF_mus.csv" \
  -t "$COUNT_OUT_DIR/transcript_count_matrix_GTF_mus.csv"
