import argparse
import csv
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from scipy.stats import fisher_exact, chi2_contingency
import pandas as pd

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

def load_dcc(dcc_path, genome_col, dcc_col):
    suffix = Path(dcc_path).suffix.lower()
    if suffix in (".xlsx", ".xls"):
        df = pd.read_excel(dcc_path, dtype=str)
    elif suffix == ".csv":
        df = pd.read_csv(dcc_path, dtype=str)
    else:
        df = pd.read_csv(dcc_path, sep="\t", dtype=str)

    if genome_col not in df.columns or dcc_col not in df.columns:
        raise SystemExit(
            f"Columns {genome_col!r}/{dcc_col!r} not found. "
            f"Available columns: {list(df.columns)}"
        )

    genome_dcc = {}
    for _, row in df.iterrows():
        gid = row[genome_col]
        dcc = row[dcc_col]
        if pd.isna(gid) or pd.isna(dcc):
            continue
        genome_dcc[str(gid).strip()] = str(dcc).strip()
    return genome_dcc

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
    ap.add_argument("--dcc-tsv", required=True, help="DCC assignment file (.xlsx, .csv, or .tsv)")
    ap.add_argument("--dcc-genome-col", default="genome_id")
    ap.add_argument("--dcc-col", default="dcc")
    ap.add_argument("--min-carriers", type=int, default=3)
    ap.add_argument("--output", required=True)
    ap.add_argument("--matrix-output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)
    genome_types = load_genome_types(args.plasmid_types_tsv, genome_ids)
    genome_dcc = load_dcc(args.dcc_tsv, args.dcc_genome_col, args.dcc_col)

    matched_genomes = [g for g in genome_types if g in genome_dcc]
    all_dcc_values = sorted(set(genome_dcc[g] for g in matched_genomes))
    all_types = sorted(set(t for types in genome_types.values() for t in types))

    print(f"{len(matched_genomes)} genomes have both a plasmid type and a DCC assignment")
    print(f"{len(all_dcc_values)} distinct DCC groups, {len(all_types)} distinct plasmid types")
    print()

    results = []
    for plasmid_type in all_types:
        carriers = [g for g in matched_genomes if plasmid_type in genome_types[g]]
        if len(carriers) < args.min_carriers:
            continue

        contingency = defaultdict(lambda: [0, 0])
        for g in matched_genomes:
            dcc = genome_dcc[g]
            has_type = plasmid_type in genome_types[g]
            contingency[dcc][0 if has_type else 1] += 1

        table = [[contingency[dcc][0], contingency[dcc][1]] for dcc in all_dcc_values]

        try:
            chi2, p_chi2, dof, expected = chi2_contingency(table)
        except ValueError:
            p_chi2 = None

        best_dcc, best_p_fisher, best_frac_in_dcc, best_direction, best_odds_ratio = None, None, None, None, None
        for dcc in all_dcc_values:
            in_dcc_with_type = contingency[dcc][0]
            in_dcc_without_type = contingency[dcc][1]
            outside_with_type = len(carriers) - in_dcc_with_type
            outside_without_type = (len(matched_genomes) - len(carriers)) - in_dcc_without_type
            two_by_two = [[in_dcc_with_type, in_dcc_without_type],
                          [outside_with_type, outside_without_type]]
            odds_ratio, p_fisher = fisher_exact(two_by_two)
            frac_in_dcc = in_dcc_with_type / len(carriers) if carriers else 0
            n_in_dcc_total = in_dcc_with_type + in_dcc_without_type
            frac_dcc_that_carry = in_dcc_with_type / n_in_dcc_total if n_in_dcc_total else 0
            direction = "enriched" if frac_dcc_that_carry > (len(carriers) / len(matched_genomes)) else "depleted"
            if best_p_fisher is None or p_fisher < best_p_fisher:
                best_dcc, best_p_fisher, best_frac_in_dcc = dcc, p_fisher, frac_in_dcc
                best_direction, best_odds_ratio = direction, odds_ratio

        results.append({
            "plasmid_type": plasmid_type,
            "n_carriers": len(carriers),
            "chi2_p": p_chi2,
            "most_associated_dcc": best_dcc,
            "direction": best_direction,
            "odds_ratio": best_odds_ratio,
            "fisher_p_vs_that_dcc": best_p_fisher,
            "frac_carriers_in_that_dcc": best_frac_in_dcc,
        })

    fisher_ps = [r["fisher_p_vs_that_dcc"] for r in results]
    adjusted = benjamini_hochberg(fisher_ps)
    for r, adj_p in zip(results, adjusted):
        r["fisher_p_fdr"] = adj_p

    results.sort(key=lambda r: r["fisher_p_fdr"])

    with open(args.output, "w") as out_fh:
        out_fh.write("plasmid_type\tn_carriers\tchi2_p\tmost_associated_dcc\tdirection\todds_ratio\tfisher_p_vs_that_dcc\tfisher_p_fdr\tfrac_carriers_in_that_dcc\n")
        for r in results:
            out_fh.write(
                f"{r['plasmid_type']}\t{r['n_carriers']}\t{r['chi2_p']}\t"
                f"{r['most_associated_dcc']}\t{r['direction']}\t{r['odds_ratio']}\t{r['fisher_p_vs_that_dcc']}\t"
                f"{r['fisher_p_fdr']}\t{r['frac_carriers_in_that_dcc']:.3f}\n"
            )

    n_significant = sum(1 for r in results if r["fisher_p_fdr"] < 0.05)
    print(f"{n_significant} plasmid types significantly associated with a specific DCC (FDR < 0.05)")
    print(f"Output written to {args.output}")

    with open(args.matrix_output, "w") as out_fh:
        out_fh.write("plasmid_type\tdcc\tn_carriers_in_dcc\tn_total_in_dcc\tfrac_carriers_of_type_in_dcc\n")
        for plasmid_type in all_types:
            carriers = [g for g in matched_genomes if plasmid_type in genome_types[g]]
            if len(carriers) < args.min_carriers:
                continue
            for dcc in all_dcc_values:
                genomes_in_dcc = [g for g in matched_genomes if genome_dcc[g] == dcc]
                carriers_in_dcc = [g for g in genomes_in_dcc if plasmid_type in genome_types[g]]
                frac = len(carriers_in_dcc) / len(carriers) if carriers else 0
                out_fh.write(f"{plasmid_type}\t{dcc}\t{len(carriers_in_dcc)}\t{len(genomes_in_dcc)}\t{frac:.4f}\n")
    print(f"Full type x DCC matrix written to {args.matrix_output}")

if __name__ == "__main__":
    main()
