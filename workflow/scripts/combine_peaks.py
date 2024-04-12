import pandas as pd

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--rev", type=str)
argparser.add_argument("--fwd", type=str)
argparser.add_argument("--sample-name", type=str)
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

fwd = pd.read_csv(
    args.fwd,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "counts", "fraction"],
)
fwd["side"] = "forward"

rev = pd.read_csv(
    args.rev,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "counts", "fraction"],
)
rev["side"] = "reverse"

merged = pd.concat([fwd, rev]).sort_values(["chrom", "start", "end"])
merged["sample"] = args.sample_name

merged.to_csv(args.output, sep="\t", header=False, index=False)
