#!/usr/bin/env bash
set -euo pipefail

THREADS=8

GTF="data/reference/gencode.vM36.annotation.gtf"

mkdir -p results/featureCounts

echo "=== STAR ==="

featureCounts \
    -T $THREADS \
    -p \
    --countReadPairs \
    -B \
    -C \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o results/featureCounts/star_counts.txt \
    results/alignments/star/*Aligned.sortedByCoord.out.bam

echo "=== HISAT2 ==="

featureCounts \
    -T $THREADS \
    -p \
    --countReadPairs \
    -B \
    -C \
    -a $GTF \
    -t exon \
    -g gene_id \
    -o results/featureCounts/hisat2_counts.txt \
    results/alignments/hisat2/*.bam

echo "=== DONE ==="
