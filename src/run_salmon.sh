#!/usr/bin/env bash
set -euo pipefail

# =========================
# CONFIGURACIÓN
# =========================

INPUT_DIR="results/nextflow/fastp"
INDEX="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.salmon"
OUTDIR="results/salmon"
LOGDIR="logs/salmon"
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

    /usr/bin/time -v salmon quant \
        -i $INDEX \
        -l A \
        -1 $R1 \
        -2 $R2 \
        -p $THREADS \
        --validateMappings \
        --gcBias \
        -o ${OUTDIR}/${base} \
        2> ${LOGDIR}/${base}.log

    echo "Finalizado $base"

done

echo "Salmon terminado"