#!/usr/bin/env bash
set -euo pipefail

# =========================
# CONFIGURACIÓN
# =========================

INPUT_DIR="results/nextflow/fastp"
GENOME_DIR="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.star"
OUTDIR="results/alignments/star"
LOGDIR="logs/star"
THREADS=8

# =========================
# DIRECTORIOS
# =========================

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

# =========================
# LOOP PRINCIPAL
# =========================

for R1 in $INPUT_DIR/*_1.clean.fastq.gz
do

    base=$(basename "$R1" _1.clean.fastq.gz)

    R2="$INPUT_DIR/${base}_2.clean.fastq.gz"

    echo "Procesando $base"

    /usr/bin/time -v STAR \
        --runThreadN $THREADS \
        --genomeDir $GENOME_DIR \
        --readFilesIn $R1 $R2 \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix ${OUTDIR}/${base}_ \
        > ${LOGDIR}/${base}.log 2>&1

    echo "Finalizado $base"

done

echo "STAR terminado"