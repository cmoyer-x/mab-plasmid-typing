import argparse
import csv
from collections import defaultdict
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
    counts = defaultdict(int)
    for r in SeqIO.parse(all_genes_faa, "fasta"):
        unit_id = r.id.rsplit("|gene", 1)[0]
        counts[unit_id] += 1
    return counts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final-labels-tsv", required=True)
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--all-genes-faa", required=True)
    ap.add_argument("--similarity-matrix", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)
    lengths = load_unit_lengths(args.plasmid_fastas_dir, genome_ids)
    gene_counts = load_gene_counts(args.all_genes_faa)
    matrix = load_matrix(args.similarity_matrix)

    unclassified = []
    with open(args.final_labels_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["method"] in ("below_floor", "no_anchor_data"):
                unclassified.append(row["unit_id"])

    rows = []
    for unit_id in unclassified:
        length = lengths.get(unit_id, "NA")
        n_genes = gene_counts.get(unit_id, 0)

        best_other = None
        best_other_score = -1
        scores_to_other_unclassified = []
        for other in unclassified:
            if other == unit_id:
                continue
            peq = get_peq(matrix, unit_id, other)
            if peq is not None:
                scores_to_other_unclassified.append(peq)
                if peq > best_other_score:
                    best_other_score = peq
                    best_other = other

        avg_to_unclassified = (
            sum(scores_to_other_unclassified) / len(scores_to_other_unclassified)
            if scores_to_other_unclassified else 0
        )

        rows.append((
            unit_id, length, n_genes,
            f"{avg_to_unclassified:.4f}",
            best_other or "NA",
            f"{best_other_score:.4f}" if best_other_score >= 0 else "NA",
        ))

    rows.sort(key=lambda r: r[1] if isinstance(r[1], int) else 0)

    with open(args.output, "w") as out_fh:
        out_fh.write("unit_id\tlength_bp\tn_genes\tavg_peq_to_other_unclassified\tmost_similar_unclassified_unit\tmost_similar_unclassified_peq\n")
        for row in rows:
            out_fh.write("\t".join(str(x) for x in row) + "\n")

    lengths_found = [r[1] for r in rows if isinstance(r[1], int)]
    gene_counts_found = [r[2] for r in rows]

    print(f"{len(rows)} unclassified units characterized")
    if lengths_found:
        print(f"Length range: {min(lengths_found)}-{max(lengths_found)} bp (median {sorted(lengths_found)[len(lengths_found)//2]})")
    print(f"Gene count range: {min(gene_counts_found)}-{max(gene_counts_found)} genes")
    n_high_mutual = sum(1 for r in rows if float(r[3]) >= 0.5)
    print(f"{n_high_mutual} units have avg PEQ >= 0.5 to other unclassified units (candidates for a novel shared cluster)")
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
