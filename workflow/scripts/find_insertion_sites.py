import numpy as np
import pandas as pd
import bioframe
from Bio import SeqIO

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--peaks", type=str)
argparser.add_argument("--genome", "-g", type=str)
argparser.add_argument("--insertion-seq", type=str, default="TA")
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()
ins_seq = args.insertion_seq
peaks = (
    pd.read_csv(
        args.peaks,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "sample_name", "fraction", "side"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "fraction": float,
            "side": str,
            "sample_name": str,
        },
    )
    .sort_values(["sample_name", "chrom", "start", "end"])
    .reset_index(drop=True)
)

peaks = bioframe.cluster(
    peaks,
    min_dist=5,
    on=["sample_name"],
    return_cluster_ids=True,
)


def determine_direction(series):
    if series.shape[0] != 2:
        return "."
    if np.all(series.to_numpy() == np.asarray(["+", "-"])):
        return "-"
    elif np.all(series.to_numpy() == np.asarray(["-", "+"])):
        return "+"
    else:
        return "."


peaks["strand"] = peaks.groupby(["cluster"])["side"].transform(determine_direction)


peaks["start"] = peaks["cluster_start"]
peaks["end"] = peaks["cluster_end"]

peaks = peaks[["chrom", "start", "end", "sample_name", "fraction", "strand"]]
peaks = (
    peaks.groupby(["chrom", "start", "end", "sample_name", "strand"])["fraction"]
    .mean()
    .reset_index()
)
peaks["score"] = (peaks["fraction"] * 1000).round().astype(int)
peaks = peaks.drop_duplicates().reset_index(drop=True)
peaks = bioframe.expand(peaks, len(args.insertion_seq))

pinpointed = []
for record in SeqIO.parse(args.genome, "fasta"):
    selected = peaks[peaks["chrom"] == record.id]
    if len(peaks) == 0:
        continue
    for i, row in selected.iterrows():
        seq = record.seq[row["start"] : row["end"]].upper()
        # Just finds the first TA
        i = seq.find(ins_seq)
        if i == -1:
            row[f"{ins_seq}_found"] = False
        else:
            row[f"{ins_seq}_found"] = True
            row["start"] = row["start"] + i
            row["end"] = row["start"] + 1
        pinpointed.append(row)
pinpointed = pd.DataFrame(pinpointed)
pinpointed = pinpointed[
    ["chrom", "start", "end", "sample_name", "score", "strand", f"{ins_seq}_found"]
]

pinpointed.sort_values(["sample_name", "chrom", "start", "end"]).to_csv(
    args.output, sep="\t", index=False, header=True
)
