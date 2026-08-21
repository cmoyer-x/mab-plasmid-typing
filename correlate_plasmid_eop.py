import argparse
import csv
from collections import defaultdict
import pandas as pd
from scipy.stats import mannwhitneyu

def load_genome_ids(summary_tsv):
    genome_ids = set()
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_ids.add(row["genome_id"])
    return genome_ids

def resolve_genome_id(unit_id, genome_ids):
    if "__" in unit_id:
        candidate = unit_id.split("__", 1)[0]
        if candidate in genome_ids:
            return candidate
    best = None
    for gid in genome_ids:
        if unit_id.startswith(gid):
            if best is None or len(gid) > len(best):
                best = gid
    return best

def load_genome_types(plasmid_types_tsv, genome_ids):
    genome_types = defaultdict(set)
    with open(plasmid_types_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["source"] == "excluded":
                continue
            gid = resolve_genome_id(row["unit_id"], genome_ids)
            if gid:
                genome_types[gid].add(row["final_type_label"])
    return genome_types

def benjamini_hochberg(pvals):
    n = len(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    adjusted = [0] * n
    prev = 1.0
    for rank, (idx, p) in reversed(list(enumerate(indexed, start=1))):
        val = min(prev, p * n / rank)
        adjusted[idx] = val
        prev = val
    return adjusted

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plasmid-types-tsv", required=True)
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--cohort-master-csv", required=True)
    ap.add_argument("--strain-col", default="strain")
    ap.add_argument("--eop-col", default="eop1_fraction")
    ap.add_argument("--min-carriers", type=int, default=3)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)
    genome_types = load_genome_types(args.plasmid_types_tsv, genome_ids)

    cohort = pd.read_csv(args.cohort_master_csv, dtype=str)
    if args.strain_col not in cohort.columns or args.eop_col not in cohort.columns:
        raise SystemExit(
            f"Columns {args.strain_col!r}/{args.eop_col!r} not found. "
            f"Available columns: {list(cohort.columns)}"
        )

    genome_eop = {}
    for _, row in cohort.iterrows():
        strain = row[args.strain_col]
        eop_raw = row[args.eop_col]
        if pd.isna(strain) or pd.isna(eop_raw):
            continue
        try:
            genome_eop[str(strain).strip()] = float(eop_raw)
        except ValueError:
            continue

    matched_genomes = [g for g in genome_types if g in genome_eop]
    all_types = sorted(set(t for types in genome_types.values() for t in types))

    print(f"{len(matched_genomes)} genomes have both a plasmid type and an EOP value")
    print(f"{len(all_types)} distinct plasmid types under consideration")
    print()

    results = []
    for plasmid_type in all_types:
        carriers = [g for g in matched_genomes if plasmid_type in genome_types[g]]
        non_carriers = [g for g in matched_genomes if plasmid_type not in genome_types[g]]
        if len(carriers) < args.min_carriers or len(non_carriers) < args.min_carriers:
            continue

        carrier_eop = [genome_eop[g] for g in carriers]
        non_carrier_eop = [genome_eop[g] for g in non_carriers]

        stat, p = mannwhitneyu(carrier_eop, non_carrier_eop, alternative="two-sided")

        median_carrier = sorted(carrier_eop)[len(carrier_eop) // 2]
        median_non_carrier = sorted(non_carrier_eop)[len(non_carrier_eop) // 2]
        direction = "higher_EOP_in_carriers" if median_carrier > median_non_carrier else "lower_EOP_in_carriers"

        results.append({
            "plasmid_type": plasmid_type,
            "n_carriers": len(carriers),
            "n_non_carriers": len(non_carriers),
            "median_eop_carriers": median_carrier,
            "median_eop_non_carriers": median_non_carrier,
            "mean_eop_carriers": sum(carrier_eop) / len(carrier_eop),
            "mean_eop_non_carriers": sum(non_carrier_eop) / len(non_carrier_eop),
            "direction": direction,
            "mannwhitney_p": p,
        })

    pvals = [r["mannwhitney_p"] for r in results]
    adjusted = benjamini_hochberg(pvals)
    for r, adj_p in zip(results, adjusted):
        r["mannwhitney_p_fdr"] = adj_p

    results.sort(key=lambda r: r["mannwhitney_p_fdr"])

    with open(args.output, "w") as out_fh:
        out_fh.write(
            "plasmid_type\tn_carriers\tn_non_carriers\tmedian_eop_carriers\t"
            "median_eop_non_carriers\tmean_eop_carriers\tmean_eop_non_carriers\t"
            "direction\tmannwhitney_p\tmannwhitney_p_fdr\n"
        )
        for r in results:
            out_fh.write(
                f"{r['plasmid_type']}\t{r['n_carriers']}\t{r['n_non_carriers']}\t"
                f"{r['median_eop_carriers']:.4f}\t{r['median_eop_non_carriers']:.4f}\t"
                f"{r['mean_eop_carriers']:.4f}\t{r['mean_eop_non_carriers']:.4f}\t"
                f"{r['direction']}\t{r['mannwhitney_p']}\t{r['mannwhitney_p_fdr']}\n"
            )

    n_significant = sum(1 for r in results if r["mannwhitney_p_fdr"] < 0.05)
    print(f"{n_significant} plasmid types significantly associated with EOP (FDR < 0.05)")
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
