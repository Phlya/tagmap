import numpy as np
import pandas as pd
import bioframe
from Bio import SeqIO

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--peaks", type=str)
argparser.add_argument("--sample-name", type=str, default=".")
argparser.add_argument("--genome", "-g", type=str)
argparser.add_argument("--insertion-seq", type=str, default="TA")
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

peaks = pd.read_csv(
    args.peaks,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "sample_name", "count", "fraction", "side"],
    dtype={
        "chrom": str,
        "start": int,
        "end": int,
        "count": int,
        "fraction": float,
        "side": str,
        "sample_name": str,
    },
)

f = peaks[peaks["side"] == "+"]
r = peaks[peaks["side"] == "-"]


overlap = bioframe.overlap(
    bioframe.expand(f, 500), bioframe.expand(r, 500), how="inner"
)
if len(overlap) == 0:
    open(args.output, "w").close()
    exit()

overlap["count"] = overlap["count"] + overlap["count_"]
overlap = overlap.drop(columns=["count_"])

overlap["score"] = (overlap["fraction"] + overlap["fraction_"]) / 2 * 1000
overlap = overlap.drop(columns=["fraction_"])

start = overlap[["start", "start_"]].max(axis=1)
end = overlap[["end", "end_"]].min(axis=1)

strand = np.where(start == overlap["start"], "+", "-")
overlap = overlap.assign(start=start, end=end, strand=strand)
overlap = overlap.drop(columns=["start_", "end_"])
# overlap = bioframe.expand(overlap, -500)

pinpointed = []
for record in SeqIO.parse(args.genome, "fasta"):
    overlaps = overlap[overlap["chrom"] == record.id]
    if len(overlaps) == 0:
        continue
    for i, row in overlaps.iterrows():
        seq = record.seq[row["start"] : row["end"]].upper()
        # Just finds the first TA
        i = seq.find(args.insertion_seq)
        if i == -1:
            row["pinpointed"] = False
        else:
            row["pinpointed"] = True
            row["start"] = row["start"] + i
            row["end"] = row["start"] + 1
        pinpointed.append(row)
pinpointed = pd.DataFrame(pinpointed)
pinpointed["name"] = args.sample_name
pinpointed = pinpointed[
    ["chrom", "start", "end", "name", "score", "strand", "pinpointed"]
]

pinpointed.sort_values(["chrom", "start", "end"]).to_csv(
    args.output, sep="\t", index=False, header=False
)
