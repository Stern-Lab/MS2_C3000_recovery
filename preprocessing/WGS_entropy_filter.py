"""
Removes indels in low-complexity sequence from a LoFreq VCF.

For each indel, Shannon entropy (bits) is computed for the reference sequence downstream (3 x indel length)
and upstream (2 x indel length) of the site; indels whose product of the two entropies is below the cutoff
are removed. SNVs are kept unchanged.

Usage:
    python WGS_entropy_filter.py --vcf in.vcf --ref reference.fasta --out final.vcf --summary summary.csv [--cutoff 1.0]
"""
import argparse
from collections import Counter

import numpy as np
import pandas as pd
import pysam


def shannon(s):
    if not s:
        return 0.0
    p = np.array(list(Counter(s).values())) / len(s)
    return float(-(p * np.log2(p)).sum())


def low_complexity(rec, seq, cutoff):
    if len(rec.ref) == len(rec.alts[0]):
        return False
    pos = rec.pos - 1
    n = max(len(rec.ref) - 1, 1)
    fwd = seq[pos + 1: pos + 1 + n * 3]
    rev = seq[max(pos - n * 2, 0): pos + 1]
    return shannon(fwd) * shannon(rev) < cutoff


parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--ref", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--summary", required=True)
parser.add_argument("--cutoff", type=float, default=1.0)
args = parser.parse_args()

ref = pysam.FastaFile(args.ref)
seqs = {name: ref.fetch(name).upper() for name in ref.references}

vcf_in = pysam.VariantFile(args.vcf)
for name, seq in seqs.items():
    if name not in vcf_in.header.contigs:
        vcf_in.header.contigs.add(name, length=len(seq))
vcf_out = pysam.VariantFile(args.out, "w", header=vcf_in.header)

stats = []
for rec in vcf_in:
    status = "FAIL_ENTROPY" if low_complexity(rec, seqs[rec.chrom], args.cutoff) else "PASS"
    if status == "PASS":
        vcf_out.write(rec)
    stats.append({"CHROM": rec.chrom, "POS": rec.pos, "REF": rec.ref, "ALT": rec.alts[0], "STATUS": status})

vcf_in.close()
vcf_out.close()
pd.DataFrame(stats).to_csv(args.summary, index=False)
