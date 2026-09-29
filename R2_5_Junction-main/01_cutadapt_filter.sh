#!/bin/bash
# Two-step cutadapt filter: keep reads with both forward primer+rDNA and reverse primer+insert.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DATA="$SCRIPT_DIR/raw_data"
OUTDIR="$SCRIPT_DIR/cutadapt_filtered"
mkdir -p "$OUTDIR"

for FASTQ in "$RAW_DATA"/*.fastq; do
    SAMPLE=$(basename "$FASTQ" .fastq)
    echo "=== Processing sample: $SAMPLE ==="

    # Step 1a: keep reads containing forward primer + rDNA sequence
    python3 -m cutadapt \
        -e 0.05 --overlap 30 --revcomp --no-indels \
        -b CGGTGGTCATGGAAGTCGAAATCCGCTAAG \
        --discard-untrimmed --action=none \
        -o "$OUTDIR/${SAMPLE}_with_F1_25S.fastq" \
        --json="$OUTDIR/${SAMPLE}_step1a.json" \
        "$FASTQ"

    # Step 1b: keep reads also containing reverse primer + insert sequence
    python3 -m cutadapt \
        -e 0.05 --overlap 30 --revcomp --no-indels \
        -b AACAAATTCCTCGGTAAAAAGCCCCaggtc \
        --discard-untrimmed --action=none \
        -o "$OUTDIR/${SAMPLE}_filtered.fastq" \
        --json="$OUTDIR/${SAMPLE}_step1b.json" \
        "$OUTDIR/${SAMPLE}_with_F1_25S.fastq"

    echo "  Output: $OUTDIR/${SAMPLE}_filtered.fastq"
done

echo ""
echo "Cutadapt filtering complete. Filtered files in: $OUTDIR/"
