"""Turn a set of Sanger .ab1 traces into one gzipped fastq.

Each read is named after its trace file, so that plate/well/direction can be
recovered from the file name later on. Also accepts fastq input, in which case
the reads are passed through unchanged apart from being concatenated.

Two basecallers are available:

* ``abi-trim`` uses the base calls embedded in the trace by the sequencer and
  applies Mott trimming, via Biopython. No external binary needed.
* ``tracy`` re-calls the bases from the trace itself, which can do better at
  ambiguous positions.
"""

import argparse
import gzip
import os
import shutil
import subprocess
import tempfile

from Bio import SeqIO

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("traces", nargs="+", help=".ab1 trace files, or fastq files")
argparser.add_argument(
    "--basecaller",
    choices=["abi-trim", "tracy"],
    default="abi-trim",
    help="Only applies to .ab1 input",
)
argparser.add_argument(
    "--tracy-args",
    default="",
    help="Extra arguments passed on to 'tracy basecall'",
)
argparser.add_argument("--output", "-o", required=True, help="Gzipped fastq")


def read_name(path):
    """Read name from the file name: 'A01_plate3.ab1' -> 'A01_plate3'."""
    name = os.path.basename(path)
    for suffix in (".gz", ".ab1", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def basecall_tracy(path, tracy_args):
    """Run tracy on one trace and return the records it called."""
    with tempfile.TemporaryDirectory() as tmpdir:
        out = os.path.join(tmpdir, "read.fastq")
        command = ["tracy", "basecall", "-f", "fastq", "-o", out]
        command += tracy_args.split()
        command += [path]
        subprocess.run(command, check=True, capture_output=True)
        # tracy appends .gz when it compresses; accept either.
        if not os.path.exists(out) and os.path.exists(out + ".gz"):
            with gzip.open(out + ".gz", "rt") as f_in, open(out, "w") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return list(SeqIO.parse(out, "fastq"))


def is_fastq(path):
    return path.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def records_for(path, basecaller, tracy_args):
    if is_fastq(path):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            return list(SeqIO.parse(f, "fastq"))
    if basecaller == "tracy":
        return basecall_tracy(path, tracy_args)
    # 'abi-trim' applies Mott's trimming algorithm to the embedded base calls.
    return list(SeqIO.parse(path, "abi-trim"))


if __name__ == "__main__":
    args = argparser.parse_args()

    n_written = 0
    with gzip.open(args.output, "wt") as out:
        for path in args.traces:
            records = records_for(path, args.basecaller, args.tracy_args)
            if not records:
                print(f"No sequence basecalled from {path}, skipping")
                continue
            name = read_name(path)
            for i, record in enumerate(records):
                # One trace is one read, so name it after the file - that is
                # where the plate, well and primer live. Reads that came from a
                # fastq keep the names they already have.
                if not is_fastq(path):
                    record.id = name if len(records) == 1 else f"{name}_{i}"
                    record.name = record.id
                record.description = ""
                SeqIO.write(record, out, "fastq")
                n_written += 1

    print(f"Wrote {n_written} reads from {len(args.traces)} files to {args.output}")
