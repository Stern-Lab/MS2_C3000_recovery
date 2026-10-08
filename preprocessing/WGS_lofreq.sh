#!/bin/bash
# LoFreq variant calling and entropy-based indel filtering for one isolate.
# Usage: WGS_lofreq.sh -i <sample.sorted.bam> -r <reference.fasta> -o <output_dir> [-a AF] [-s SB_ALPHA] [-c COV_MIN] [-e ENTROPY_CUTOFF] [-q INDEL_QUAL]
set -euo pipefail

AF=0.01
SB_ALPHA=0.01
COV_MIN=30
ENTROPY_CUTOFF=1.0
INDEL_QUAL=16
BAM_INPUT=""
REF=""
OUT_DIR=""

usage() { grep "^# Usage" "$0" | cut -c3-; exit 1; }

while getopts "i:r:o:a:s:c:e:q:h" opt; do
    case $opt in
        i) BAM_INPUT=$OPTARG ;;
        r) REF=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        a) AF=$OPTARG ;;
        s) SB_ALPHA=$OPTARG ;;
        c) COV_MIN=$OPTARG ;;
        e) ENTROPY_CUTOFF=$OPTARG ;;
        q) INDEL_QUAL=$OPTARG ;;
        *) usage ;;
    esac
done
[ -z "$BAM_INPUT" ] || [ -z "$REF" ] || [ -z "$OUT_DIR" ] && usage

mkdir -p "$OUT_DIR"
SAMPLE=$(basename "$BAM_INPUT" .sorted.bam)
OUT="$OUT_DIR/$SAMPLE"

lofreq indelqual --uniform "$INDEL_QUAL" -f "$REF" -o "${OUT}_with_quals.bam" "$BAM_INPUT"
samtools index "${OUT}_with_quals.bam"

lofreq call --no-default-filter --call-indels -f "$REF" -o "${OUT}_raw.vcf" "${OUT}_with_quals.bam"

lofreq filter -i "${OUT}_raw.vcf" -o "${OUT}_lofreq_filtered.vcf" \
    -a "$AF" --sb-alpha "$SB_ALPHA" --sb-incl-indels --cov-min "$COV_MIN"

python3 "$(dirname "$0")/WGS_entropy_filter.py" --vcf "${OUT}_lofreq_filtered.vcf" --ref "$REF" \
    --cutoff "$ENTROPY_CUTOFF" --out "${OUT}_final.vcf" --summary "${OUT}_summary.csv"
