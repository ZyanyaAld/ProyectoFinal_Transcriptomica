#!/usr/bin/env bash
set -euo pipefail

# =========================
# CONFIGURACIÓN
# =========================

INPUT_DIR="results/nextflow/fastp"
INDEX="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat"
OUTDIR="results/alignments/hisat2"
LOGDIR="logs/hisat2"
THREADS=8

# =========================
# DIRECTORIOS
# =========================

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

# =========================
# LOOP
# =========================

for R1 in $INPUT_DIR/*_1.clean.fastq.gz
do

    base=$(basename "$R1" _1.clean.fastq.gz)

    R2="$INPUT_DIR/${base}_2.clean.fastq.gz"

    echo "Procesando $base"

    /usr/bin/time -v hisat2 \
        -p $THREADS \
        -x $INDEX \
        -1 $R1 \
        -2 $R2 \
        --summary-file ${OUTDIR}/${base}_summary.txt \
        2> ${LOGDIR}/${base}.log \
        | samtools sort -@ $THREADS \
        -o ${OUTDIR}/${base}.bam

    echo "Finalizado $base"

done

echo "HISAT2 terminado"