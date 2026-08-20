import argparse
import csv
from pathlib import Path
from Bio import SeqIO

PAPER_CLUSTERS = {
    "pA": [("GD08", 9547), ("GD42", 9547)],
    "pB": [
        ("GD18", 25000), ("GD62", 25000), ("GD69", 25000), ("GD95", 25000),
        ("GD108A", 25000), ("GD108B", 25000), ("GD23", 25002), ("GD36", 24995),
        ("GD47", 24995), ("GD42", 24993), ("GD72", 24985), ("GD87", 24994),
    ],
    "pC": [
        ("GD22", 18117), ("GD24", 18117), ("GD34", 18117), ("GD75", 18117),
        ("GD100A", 18117), ("GD100B", 18117), ("GD39", 18117),
        ("GD62", 18612), ("GD69", 18611), ("GD95", 18611),
    ],
    "pD": [("GD19", 18605), ("GD45", 19406), ("GD85", 23374)],
    "pE": [("GD33", 25996), ("GD36", 24259)],
    "pF": [("GD02", 30963), ("GD25", 31413), ("GD54", 31413), ("GD102", 31413), ("GD86", 31343)],
    "pG": [("GD25", 27424), ("GD45", 27427), ("GD86", 27424), ("GD102", 27425)],
    "pH": [("GD58", 92821)],
    "singleton_pGD13": [("GD13", 21881)],
    "singleton_pGD21-1": [("GD21", 112633)],
    "singleton_pGD21-2": [("GD21", 155609)],
    "singleton_pGD22-1": [("GD22", 19694)],
    "singleton_pGD25-3": [("GD25", 23599)],
    "singleton_pGD51": [("GD51", 23656)],
    "singleton_pGD52": [("GD52", 22216)],
    "singleton_pGD104": [("GD104", 96413)],
    "singleton_pATCC19977": [],
}

def load_genome_ids(summary_tsv):
    genome_ids = set()
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_ids.add(row["genome_id"])
    return genome_ids

def load_cluster_assignments(cluster_tsv):
    rows = []
    with open(cluster_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows

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

def load_unit_lengths(plasmid_fastas_dir, genome_ids):
    lengths = {}
    plasmid_dir = Path(plasmid_fastas_dir)
    for genome_id in genome_ids:
        fasta_path = plasmid_dir / f"{genome_id}_plasmids.fasta"
        if not fasta_path.exists():
            continue
        for r in SeqIO.parse(fasta_path, "fasta"):
            unit_id = r.id if r.id.startswith(genome_id) else f"{genome_id}__{r.id}"
            lengths[unit_id] = len(r.seq)
    return lengths

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--cluster-assignments", required=True)
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--length-tolerance", type=float, default=0.10)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)
    cluster_rows = load_cluster_assignments(args.cluster_assignments)
    unit_lengths = load_unit_lengths(args.plasmid_fastas_dir, genome_ids)

    unit_to_cluster = {}
    for row in cluster_rows:
        unit_to_cluster[row["unit_id"]] = row["cluster_id"]

    genome_to_units = {}
    for unit_id in unit_to_cluster:
        gid = resolve_genome_id(unit_id, genome_ids)
        if gid:
            genome_to_units.setdefault(gid, []).append(unit_id)

    crosswalk = {}
    not_found = []
    unit_claimed_by = {}

    for paper_cluster, strain_lengths in PAPER_CLUSTERS.items():
        found_any = False
        for strain, expected_length in strain_lengths:
            candidates = genome_to_units.get(strain, [])
            if not candidates:
                continue
            best_unit = None
            best_diff = None
            for unit_id in candidates:
                actual_length = unit_lengths.get(unit_id)
                if actual_length is None:
                    continue
                diff = abs(actual_length - expected_length) / expected_length
                if diff <= args.length_tolerance and (best_diff is None or diff < best_diff):
                    best_unit = unit_id
                    best_diff = diff
            if best_unit is None:
                continue
            found_any = True
            phamclust_cluster = unit_to_cluster[best_unit]
            crosswalk.setdefault(phamclust_cluster, {}).setdefault(paper_cluster, []).append(
                (strain, best_unit, unit_lengths[best_unit], expected_length)
            )
            unit_claimed_by.setdefault(best_unit, []).append(paper_cluster)
        if not found_any and strain_lengths:
            not_found.append(paper_cluster)

    with open(args.output, "w") as out_fh:
        out_fh.write("phamclust_cluster\tpaper_cluster\tevidence_strain\tevidence_unit_id\tactual_length\texpected_length\n")
        for phamclust_cluster, paper_hits in sorted(crosswalk.items()):
            for paper_cluster, evidence in paper_hits.items():
                for strain, unit_id, actual_len, expected_len in evidence:
                    out_fh.write(f"{phamclust_cluster}\t{paper_cluster}\t{strain}\t{unit_id}\t{actual_len}\t{expected_len}\n")

    print(f"Crosswalk written to {args.output}")
    print()
    print("phamclust cluster -> paper cluster(s):")
    for phamclust_cluster, paper_hits in sorted(crosswalk.items()):
        labels = ",".join(paper_hits.keys())
        n_paper_labels = len(paper_hits)
        flag = "  <-- AMBIGUOUS: maps to multiple paper clusters" if n_paper_labels > 1 else ""
        print(f"  {phamclust_cluster}: {labels}{flag}")
    print()
    multi_claimed = {u: c for u, c in unit_claimed_by.items() if len(set(c)) > 1}
    if multi_claimed:
        print("Units claimed by more than one paper cluster (worth checking):")
        for unit_id, claims in multi_claimed.items():
            print(f"  {unit_id}: {claims}")
        print()
    if not_found:
        print(f"Paper clusters with no representative strain found in this cohort: {', '.join(not_found)}")

if __name__ == "__main__":
    main()

