from pathlib import Path
import time
import numpy as np

import bioframe
import pandas as pd
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", "-i", type=str)
argparser.add_argument("--chromsizes", type=str)
argparser.add_argument("--name", type=str, default=None)
argparser.add_argument("--output", "-o", type=str)
argparser.add_argument("--output-bigwig", type=str, default=None, required=False)


if __name__ == "__main__":
    args = argparser.parse_args()

    bg = pd.read_table(
        args.input,
        header=None,
        names=["chrom", "start", "end", "value", "coverage"],
        dtype={
            "chrom": str,
            "start": np.int64,
            "end": np.int64,
            "value": np.int64,
            "coverage": np.float32,
        },
    )
    chromsizes = bioframe.read_chromsizes(args.chromsizes)
    bg = bioframe.trim(bg, chromsizes).dropna()
    bg[["start", "end"]] = bg[["start", "end"]].astype(int)
    with open(args.output, "w") as f:
        f.write(f"track name={args.name} type=bedGraph\n")
    bg[["chrom", "start", "end", "value"]].to_csv(
        args.output, sep="\t", index=False, header=False, mode="a"
    )

    if args.output_bigwig is not None:
        bioframe.to_bigwig(
            bg,
            chromsizes,
            args.output_bigwig,
        )
