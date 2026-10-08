"""
Mean chromosome and plasmid read depth per isolate, used for Figure 4B.

Input: one sorted, indexed BAM per sample named <sample>.sorted.bam, with sample names as in PRJNA1536727
(c1-c10, i1-i20). Reads were aligned to the E. coli ATCC 15597 reference genome with
    minimap2 -ax map-ont <reference> <sample>.fastq.gz | samtools sort -o <sample>.sorted.bam
(minimap2 2.30, samtools 1.22). Depth is taken from the pysam pileup (secondary, duplicate, QC-fail and
unmapped reads excluded); the larger contig is the chromosome and the smaller the F plasmid.

Usage (from the repository root; the BAM folder is taken from `wgs_bam_dir` in config.yaml):
    python preprocessing/WGS_coverage.py
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import yaml

ROOT = Path(__file__).resolve().parents[1]
config = yaml.safe_load(open(ROOT / "config.yaml"))

parser = argparse.ArgumentParser()
parser.add_argument("--bam-dir", type=Path, default=Path(config["wgs_bam_dir"]))
parser.add_argument("--out", type=Path, default=ROOT / "data" / "WGS_plasmid_chromosome_coverage.csv")
args = parser.parse_args()

rows = []
for bam_path in sorted(args.bam_dir.glob("*.sorted.bam")):
    sample = bam_path.name[: -len(".sorted.bam")]
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contigs = sorted(bam.references, key=bam.get_reference_length, reverse=True)
        if len(contigs) != 2:
            raise ValueError(f"{sample}: expected 2 contigs, found {contigs}")
        mean_depth = {}
        for contig in contigs:
            depth = np.zeros(bam.get_reference_length(contig), dtype=np.int64)
            for col in bam.pileup(contig):
                depth[col.reference_pos] = col.nsegments
            mean_depth[contig] = depth.mean()
    rows.append({"sample": sample,
                 "chrom_coverage": mean_depth[contigs[0]],
                 "plasmid_coverage": mean_depth[contigs[1]]})
    print(f"{sample}: chromosome {mean_depth[contigs[0]]:.1f}x, plasmid {mean_depth[contigs[1]]:.1f}x")

df = pd.DataFrame(rows)
args.out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(args.out, index=False)
print(f"{len(df)} samples written to {args.out}")
