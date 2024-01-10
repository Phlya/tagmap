import pandas as pd
import bioframe

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--forward-peaks", "-f", type=str)
argparser.add_argument("--reverse-peaks", "-r", type=str)
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

f = pd.read_csv(
    args.forward_peaks,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "name", "count"],
)
r = pd.read_csv(
    args.reverse_peaks,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "name", "count"],
)

overlap = bioframe.overlap(f, r, how="inner")
overlap = overlap[["chrom", "start_", "end"]]
overlap.columns = ["chrom", "start", "end"]
overlap.to_csv(args.output, sep="\t", index=False, header=False)
