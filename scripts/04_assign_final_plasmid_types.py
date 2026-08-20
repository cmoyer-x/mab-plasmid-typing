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

PHIX174_LENGTH = 5386
PHIX174_LENGTH_TOLERANCE = 0.02

def is_likely_phix(length):
    if length is None:
        return False
    return abs(length - PHIX174_LENGTH) / PHIX174_LENGTH <= PHIX174_LENGTH_TOLERANCE

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

def load_gene_counts(all_genes_faa):
    counts = {}
    for r in SeqIO.parse(all_genes_faa, "fasta"):
        unit_id = r.id.rsplit("|gene", 1)[0]
        counts[unit_id] = counts.get(unit_id, 0) + 1
    return counts

def load_matrix(matrix_tsv):
    with open(matrix_tsv) as fh:
        lines = fh.read().strip().split("\n")
    header = lines[0].split("\t")
    col_ids = header[1:]
    matrix = {}
    for line in lines[1:]:
        fields = line.split("\t")
        row_id = fields[0]
        values = fields[1:]
        matrix[row_id] = dict(zip(col_ids, values))
    return matrix

def get_peq(matrix, a, b):
    if a in matrix and b in matrix[a]:
        return float(matrix[a][b])
    if b in matrix and a in matrix[b]:
        return float(matrix[b][a])
    return None

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

def find_anchors(unit_ids, genome_ids, lengths, length_tolerance):
    genome_to_units = {}
    for unit_id in unit_ids:
        gid = resolve_genome_id(unit_id, genome_ids)
        if gid:
            genome_to_units.setdefault(gid, []).append(unit_id)

    anchor_label = {}
    for paper_cluster, strain_lengths in PAPER_CLUSTERS.items():
        for strain, expected_length in strain_lengths:
            candidates = genome_to_units.get(strain, [])
            best_unit, best_diff = None, None
            for unit_id in candidates:
                actual_length = lengths.get(unit_id)
                if actual_length is None:
                    continue
                diff = abs(actual_length - expected_length) / expected_length
                if diff <= length_tolerance and (best_diff is None or diff < best_diff):
                    best_unit, best_diff = unit_id, diff
            if best_unit is not None:
                anchor_label.setdefault(best_unit, set()).add(paper_cluster)
    return anchor_label

def cluster_by_threshold(unit_ids, matrix, threshold):
    parent = {u: u for u in unit_ids}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, a in enumerate(unit_ids):
        for b in unit_ids[i + 1:]:
            peq = get_peq(matrix, a, b)
            if peq is not None and peq >= threshold:
                union(a, b)

    clusters = {}
    for u in unit_ids:
        clusters.setdefault(find(u), []).append(u)
    return clusters

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--all-genes-faa", required=True)
    ap.add_argument("--similarity-matrix", required=True)
    ap.add_argument("--length-tolerance", type=float, default=0.10)
    ap.add_argument("--nearest-floor", type=float, default=0.50)
    ap.add_argument("--nearest-margin", type=float, default=0.05)
    ap.add_argument("--novel-clu-thresh", type=float, default=0.50)
    ap.add_argument("--min-length", type=int, default=3000)
    ap.add_argument("--min-genes", type=int, default=4)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)
    lengths = load_unit_lengths(args.plasmid_fastas_dir, genome_ids)
    gene_counts = load_gene_counts(args.all_genes_faa)
    matrix = load_matrix(args.similarity_matrix)
    unit_ids = list(matrix.keys())

    anchor_label = find_anchors(unit_ids, genome_ids, lengths, args.length_tolerance)

    anchor_group_members = {}
    for unit_id, labels in anchor_label.items():
        for label in labels:
            anchor_group_members.setdefault(label, set()).add(unit_id)

    results = {}

    for unit_id in unit_ids:
        if unit_id in anchor_label:
            label = "+".join(sorted(anchor_label[unit_id]))
            results[unit_id] = (label, "paper_cluster_anchor", "NA")
            continue

        label_avgs = {}
        for label, members in anchor_group_members.items():
            scores = [get_peq(matrix, unit_id, m) for m in members if m != unit_id]
            scores = [s for s in scores if s is not None]
            if scores:
                label_avgs[label] = sum(scores) / len(scores)

        if not label_avgs:
            results[unit_id] = (None, "pending_novel", "NA")
            continue

        ranked = sorted(label_avgs.items(), key=lambda kv: kv[1], reverse=True)
        best_label, best_score = ranked[0]
        second_score = ranked[1][1] if len(ranked) > 1 else None

        if best_score < args.nearest_floor:
            results[unit_id] = (None, "pending_novel", f"{best_score:.4f}")
        elif second_score is not None and (best_score - second_score) < args.nearest_margin:
            results[unit_id] = (best_label, "paper_cluster_ambiguous", f"{best_score:.4f}")
        else:
            results[unit_id] = (best_label, "paper_cluster_nearest", f"{best_score:.4f}")

    pending = [u for u, (label, source, _) in results.items() if source == "pending_novel"]

    excluded = []
    excluded_phix = []
    to_cluster = []
    for u in pending:
        length = lengths.get(u, 0)
        n_genes = gene_counts.get(u, 0)
        if is_likely_phix(length):
            excluded_phix.append(u)
        elif length < args.min_length or n_genes < args.min_genes:
            excluded.append(u)
        else:
            to_cluster.append(u)

    import string
    used_letters = set("ABCDEFGH")
    available_letters = [ltr for ltr in string.ascii_uppercase if ltr not in used_letters]

    novel_clusters = cluster_by_threshold(to_cluster, matrix, args.novel_clu_thresh)
    ranked_clusters = sorted(novel_clusters.values(), key=len, reverse=True)

    letter_idx = 0
    singleton_strain_counts = {}
    novel_cluster_count = 0
    for members in ranked_clusters:
        if len(members) > 1:
            if letter_idx < len(available_letters):
                label = f"p{available_letters[letter_idx]}"
            else:
                label = f"p{available_letters[letter_idx % len(available_letters)]}{letter_idx // len(available_letters) + 1}"
            letter_idx += 1
            novel_cluster_count += 1
            for m in members:
                results[m] = (label, "novel_cluster", "NA")
        else:
            unit_id = members[0]
            strain = resolve_genome_id(unit_id, genome_ids) or unit_id
            singleton_strain_counts[strain] = singleton_strain_counts.get(strain, 0) + 1
            n = singleton_strain_counts[strain]
            suffix = f"-{n}" if n > 1 else ""
            label = f"singleton_p{strain}{suffix}"
            results[unit_id] = (label, "novel_cluster", "NA")

    for u in excluded:
        results[u] = ("excluded_short_fragment", "excluded", "NA")

    for u in excluded_phix:
        results[u] = ("excluded_phix_contamination", "excluded", "NA")

    with open(args.output, "w") as out_fh:
        out_fh.write("unit_id\tfinal_type_label\tsource\tbest_avg_peq\tlength_bp\tn_genes\n")
        for unit_id in sorted(results.keys()):
            label, source, score = results[unit_id]
            length = lengths.get(unit_id, "NA")
            n_genes = gene_counts.get(unit_id, "NA")
            out_fh.write(f"{unit_id}\t{label}\t{source}\t{score}\t{length}\t{n_genes}\n")

    from collections import Counter
    source_counts = Counter(source for _, source, _ in results.values())
    print(f"Total units: {len(results)}")
    for source, count in source_counts.items():
        print(f"  {source}: {count}")
    print(f"Novel lettered clusters found: {novel_cluster_count} (pI through p{available_letters[novel_cluster_count - 1] if 0 < novel_cluster_count <= len(available_letters) else '?'})")
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
