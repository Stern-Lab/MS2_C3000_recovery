"""
Builds the MS2 strand-specific depth table and mapped-read counts used for Figure S2.

Input directory: per-sample samtools depth files (*.sorted.ms2.forward.txt / *.sorted.ms2.reverse.txt)
and filtered alignments (*.HQ_mapped.sorted.bam), as produced by the bulk RNA-seq pipeline.

Usage:
    python bulk_MS2_coverage.py --input-dir /path/to/alignments --metadata ../data/bulk_RNAseq_sample_metadata.csv --out-dir ../data
"""
import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input-dir", required=True, type=Path)
parser.add_argument("--metadata", required=True, type=Path)
parser.add_argument("--out-dir", required=True, type=Path)
parser.add_argument("--samtools", default="samtools")
args = parser.parse_args()

meta = pd.read_csv(args.metadata)
meta = meta[~meta["sample"].str.upper().str.endswith("-ON")]
samples = set(meta["sample"])
infected = set(meta.loc[meta["infected"] == "yes", "sample"])
strand_re = re.compile(r"\.(forward|reverse)\.")

tables = []
for f in sorted(args.input_dir.glob("*/*.sorted.ms2.*.txt")):
    sample = f.name.split("_")[0]
    m = strand_re.search(f.name)
    if sample not in samples or not m:
        continue
    d = pd.read_csv(f, sep="\t", names=["chromosome", "position", "coverage"])
    d["sample"] = sample
    d["strand"] = m.group(1)
    tables.append(d[["sample", "strand", "position", "coverage"]])
depth = pd.concat(tables, ignore_index=True)

reads = []
for bam in sorted(args.input_dir.glob("*/*.HQ_mapped.sorted.bam")):
    sample = bam.name.split("_")[0]
    if sample not in samples:
        continue
    n = subprocess.run([args.samtools, "view", "-c", str(bam)], capture_output=True, text=True, check=True)
    reads.append({"sample": sample, "mapped_reads": int(n.stdout.strip())})
reads = pd.DataFrame(reads)

for name, found in (("depth files", set(depth["sample"])), ("BAM files", set(reads["sample"]))):
    missing = samples - found
    if missing & infected:
        raise ValueError(f"Infected samples without {name}: {sorted(missing & infected)}")
    if missing:
        print(f"Control samples without {name} (not used in Figure S2): {sorted(missing)}")

args.out_dir.mkdir(parents=True, exist_ok=True)
depth.to_csv(args.out_dir / "bulk_MS2_strand_depth.csv.gz", index=False)
reads.to_csv(args.out_dir / "bulk_mapped_reads.csv", index=False)
print(f"{depth['sample'].nunique()} samples, {len(depth)} depth rows")
