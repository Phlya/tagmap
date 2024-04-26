import argparse
import pyfastx
import json

# Get the positions of the forward and reverse primers in the reference genome
# and save them to a yaml file

parser = argparse.ArgumentParser()
parser.add_argument("--genome", type=str)
parser.add_argument("--cassette-name", type=str)
parser.add_argument("--forward-primer", type=str)
parser.add_argument("--reverse-primer", type=str)
parser.add_argument("--output", "-o", type=str)
args = parser.parse_args()

forward_primer_sequence = args.forward_primer.upper()
reverse_primer_sequence = args.reverse_primer.upper()
complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
reverse_primer_revcomp = "".join(
    complement.get(base, base) for base in reversed(reverse_primer_sequence)
)
for name, seq in pyfastx.Fasta(args.genome, build_index=False):
    if name == "SB_cargo_CTCF":
        f_start = seq.index(forward_primer_sequence) + 1
        r_start = seq.find(reverse_primer_revcomp)
        r_end = r_start + len(reverse_primer_sequence)
        break
to_save = {"forward_ITR_primer_position": f_start, "reverse_ITR_primer_position": r_end}
with open(args.output, "w") as f:
    f.write(json.dumps(to_save))
