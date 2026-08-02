"""Summarise NGS mobilization rates and insertion-site sidedness.

Turns the per-sample pairtools stats.yml files into one table of how many
read pairs from each library actually captured a transposition event -
anchored at an ITR primer, rather than just mapping anywhere - and turns the
combined insertion sites file into a second table of how many of those sites
were seen from both sides of the cassette versus only one.
"""

import argparse
import os

import pandas as pd
import yaml

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument("--stats-yml", nargs="*", default=[])
argparser.add_argument("--sites", required=True, help="all_sites.bed")
argparser.add_argument("--output-mapping", "-o", required=True)
argparser.add_argument("--output-sidedness", required=True)
args = argparser.parse_args()

MAPPING_COLUMNS = [
    "sample_name",
    "total_pairs",
    "mapped_pairs",
    "frac_mapped",
    "forward_junction_pairs",
    "reverse_junction_pairs",
    "mobilized_pairs",
    "frac_mobilized",
]

SIDEDNESS_COLUMNS = [
    "sample_name",
    "n_sites",
    "n_two_sided",
    "frac_two_sided",
    "n_one_sided",
    "n_forward_only",
    "n_reverse_only",
]


def sample_name_from_path(path):
    return os.path.basename(path).rsplit("_stats.yml", 1)[0]


rows = []
for path in args.stats_yml:
    with open(path) as f:
        stat = yaml.safe_load(f) or {}
    total = stat.get("no_filter", {}).get("total", 0)
    mapped = stat.get("no_filter", {}).get("total_mapped", 0)
    forward = stat.get("forward", {}).get("total", 0)
    reverse = stat.get("reverse", {}).get("total", 0)
    mobilized = forward + reverse
    rows.append(
        {
            "sample_name": sample_name_from_path(path),
            "total_pairs": total,
            "mapped_pairs": mapped,
            "frac_mapped": mapped / total if total else float("nan"),
            "forward_junction_pairs": forward,
            "reverse_junction_pairs": reverse,
            "mobilized_pairs": mobilized,
            "frac_mobilized": mobilized / total if total else float("nan"),
        }
    )

mapping = pd.DataFrame(rows, columns=MAPPING_COLUMNS).sort_values("sample_name")
mapping.to_csv(args.output_mapping, sep="\t", index=False)

sites = pd.read_csv(args.sites, sep="\t", dtype={"chrom": str})

if sites.shape[0] == 0:
    sidedness = pd.DataFrame(columns=SIDEDNESS_COLUMNS)
else:

    def summarize(group):
        n = group.shape[0]
        forward_only = (group["site_sides"] == "forward_only").sum()
        reverse_only = (group["site_sides"] == "reverse_only").sum()
        two_sided = (group["site_sides"] == "both").sum()
        return pd.Series(
            {
                "n_sites": n,
                "n_two_sided": two_sided,
                "frac_two_sided": two_sided / n if n else float("nan"),
                "n_one_sided": forward_only + reverse_only,
                "n_forward_only": forward_only,
                "n_reverse_only": reverse_only,
            }
        )

    sidedness = (
        sites.groupby("sample_name")
        .apply(summarize, include_groups=False)
        .reset_index()[SIDEDNESS_COLUMNS]
    )
    # groupby.apply returns one Series per group; mixing ints and the float
    # frac_two_sided in that Series forces the whole thing to float64.
    int_columns = ["n_sites", "n_two_sided", "n_one_sided", "n_forward_only", "n_reverse_only"]
    sidedness[int_columns] = sidedness[int_columns].astype(int)

sidedness.to_csv(args.output_sidedness, sep="\t", index=False)

print(
    f"{mapping.shape[0]} NGS samples: "
    f"{int(mapping['mobilized_pairs'].sum()) if mapping.shape[0] else 0} mobilized pairs; "
    f"{sidedness['n_sites'].sum() if sidedness.shape[0] else 0} insertion sites, "
    f"{sidedness['n_two_sided'].sum() if sidedness.shape[0] else 0} seen from both sides"
)
