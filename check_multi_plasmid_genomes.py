import argparse
import csv
from collections import defaultdict

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plasmid-types-tsv", required=True)
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)

    genome_types = defaultdict(set)
    genome_units = defaultdict(list)
    with open(args.plasmid_types_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["source"] == "excluded":
                continue
            gid = resolve_genome_id(row["unit_id"], genome_ids)
            if gid:
                genome_types[gid].add(row["final_type_label"])
                genome_units[gid].append(row["unit_id"])

    multi_plasmid_genomes = {g: types for g, types in genome_types.items() if len(types) > 1}

    with open(args.output, "w") as out_fh:
        out_fh.write("genome_id\tn_distinct_types\tplasmid_types\tunit_ids\n")
        for g in sorted(genome_types, key=lambda x: -len(genome_types[x])):
            types = genome_types[g]
            out_fh.write(f"{g}\t{len(types)}\t{','.join(sorted(types))}\t{','.join(genome_units[g])}\n")

    print(f"{len(genome_types)} genomes have at least one typed plasmid")
    print(f"{len(multi_plasmid_genomes)} genomes carry more than one distinct plasmid type")
    print()
    if multi_plasmid_genomes:
        print("Genomes with multiple distinct plasmid types:")
        for g, types in sorted(multi_plasmid_genomes.items(), key=lambda x: -len(x[1])):
            print(f"  {g}: {sorted(types)}")
    print(f"Full breakdown written to {args.output}")

if __name__ == "__main__":
    main()
