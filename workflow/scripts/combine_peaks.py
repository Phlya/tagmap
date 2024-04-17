import pandas as pd
import bioframe
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--rev", type=str)
argparser.add_argument("--fwd", type=str)
argparser.add_argument("--blacklist", type=str, default=None)
argparser.add_argument("--sample-name", type=str)
argparser.add_argument("--output", "-o", type=str)
args = argparser.parse_args()

fwd = pd.read_csv(
    args.fwd,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "counts", "fraction"],
)
fwd["side"] = "+"

rev = pd.read_csv(
    args.rev,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "counts", "fraction"],
)
rev["side"] = "-"
merged = pd.concat([fwd, rev]).sort_values(["chrom", "start", "end"])

if args.blacklist is not None:
    blacklist = pd.read_csv(
        args.blacklist,
        sep="\t",
        header=None,
        comment="#",
        names=["chrom", "start", "end"],
    )
    merged = bioframe.setdiff(merged, blacklist)

merged["sample"] = args.sample_name
merged = (
    merged[["chrom", "start", "end", "sample", "counts", "fraction", "side"]]
    .sort_values(["chrom", "start", "end"])
    .to_csv(args.output, sep="\t", header=False, index=False)
)
