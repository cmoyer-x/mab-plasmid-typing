import argparse
import csv
from pathlib import Path
from Bio import SeqIO

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

def load_genome_ids(summary_tsv):
    genome_ids = set()
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_ids.add(row["genome_id"])
    return genome_ids

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
        root = find(u)
        clusters.setdefault(root, []).append(u)
    return clusters

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final-labels-tsv", required=True)
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--all-genes-faa", required=True)
    ap.add_argument("--similarity-matrix", required=True)
    ap.add_argument("--min-length", type=int, default=3000)
    ap.add_argument("--min-genes", type=int, default=4)
    ap.add_argument("--clu-thresh", type=float, default=0.5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)
    lengths = load_unit_lengths(args.plasmid_fastas_dir, genome_ids)
    gene_counts = load_gene_counts(args.all_genes_faa)
    matrix = load_matrix(args.similarity_matrix)

    final_rows = {}
    with open(args.final_labels_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            final_rows[row["unit_id"]] = row

    unclassified = [
        u for u, row in final_rows.items()
        if row["method"] in ("below_floor", "no_anchor_data")
    ]

    excluded = []
    to_cluster = []
    for u in unclassified:
        length = lengths.get(u, 0)
        n_genes = gene_counts.get(u, 0)
        if length < args.min_length or n_genes < args.min_genes:
            excluded.append(u)
        else:
            to_cluster.append(u)

    novel_clusters = cluster_by_threshold(to_cluster, matrix, args.clu_thresh)

    ranked_clusters = sorted(novel_clusters.values(), key=len, reverse=True)
    unit_to_novel_label = {}
    novel_counter = 1
    for members in ranked_clusters:
        if len(members) > 1:
            label = f"novel_{novel_counter}"
            novel_counter += 1
        else:
            label = "novel_singleton"
        for m in members:
            unit_to_novel_label[m] = (label, len(members))

    output_rows = []
    for unit_id, row in final_rows.items():
        length = lengths.get(unit_id, "NA")
        n_genes = gene_counts.get(unit_id, "NA")
        if row["method"] in ("anchor", "nearest_anchor", "ambiguous_nearest"):
            output_rows.append((unit_id, row["assigned_label"], "paper_cluster", length, n_genes))
        elif unit_id in excluded:
            output_rows.append((unit_id, "excluded_short_fragment", "excluded", length, n_genes))
        elif unit_id in unit_to_novel_label:
            label, size = unit_to_novel_label[unit_id]
            output_rows.append((unit_id, f"{label}_n{size}", "novel_cluster", length, n_genes))
        else:
            output_rows.append((unit_id, "unresolved", "unresolved", length, n_genes))

    with open(args.output, "w") as out_fh:
        out_fh.write("unit_id\tfinal_type_label\tsource\tlength_bp\tn_genes\n")
        for row in sorted(output_rows):
            out_fh.write("\t".join(str(x) for x in row) + "\n")

    n_paper = sum(1 for r in output_rows if r[2] == "paper_cluster")
    n_novel_multi = sum(1 for r in output_rows if r[2] == "novel_cluster" and "singleton" not in r[1])
    n_novel_singleton = sum(1 for r in output_rows if r[2] == "novel_cluster" and "singleton" in r[1])
    n_excluded = sum(1 for r in output_rows if r[2] == "excluded")

    print(f"Total units: {len(output_rows)}")
    print(f"  paper_cluster (anchor/nearest/ambiguous): {n_paper}")
    print(f"  novel_cluster, multi-member: {n_novel_multi}")
    print(f"  novel_cluster, singleton: {n_novel_singleton}")
    print(f"  excluded (too short/too few genes): {n_excluded}")
    print(f"Novel clusters found: {novel_counter - 1}")
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
