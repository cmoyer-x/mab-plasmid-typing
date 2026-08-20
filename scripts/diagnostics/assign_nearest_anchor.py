import argparse
import csv
from collections import defaultdict

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

def load_anchors(crosswalk_tsv):
    anchors = defaultdict(set)
    unit_to_labels = defaultdict(set)
    with open(crosswalk_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            label = row["paper_cluster"]
            unit_id = row["evidence_unit_id"]
            anchors[label].add(unit_id)
            unit_to_labels[unit_id].add(label)
    return anchors, unit_to_labels

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosswalk-tsv", required=True)
    ap.add_argument("--similarity-matrix", required=True)
    ap.add_argument("--floor", type=float, default=0.50)
    ap.add_argument("--margin", type=float, default=0.05)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = load_matrix(args.similarity_matrix)
    anchors, unit_to_labels = load_anchors(args.crosswalk_tsv)

    all_units = list(matrix.keys())

    rows = []
    for unit_id in all_units:
        if unit_id in unit_to_labels:
            labels = unit_to_labels[unit_id]
            label = "+".join(sorted(labels))
            method = "anchor" if len(labels) == 1 else "anchor_conflict"
            rows.append((unit_id, label, method, "NA", "NA", "NA"))
            continue

        label_avgs = {}
        for label, members in anchors.items():
            scores = []
            for member in members:
                if member == unit_id:
                    continue
                peq = get_peq(matrix, unit_id, member)
                if peq is not None:
                    scores.append(peq)
            if scores:
                label_avgs[label] = sum(scores) / len(scores)

        if not label_avgs:
            rows.append((unit_id, "unclassified", "no_anchor_data", "NA", "NA", "NA"))
            continue

        ranked = sorted(label_avgs.items(), key=lambda kv: kv[1], reverse=True)
        best_label, best_score = ranked[0]
        second_label, second_score = ranked[1] if len(ranked) > 1 else (None, None)

        if best_score < args.floor:
            rows.append((unit_id, "unclassified", "below_floor", f"{best_score:.4f}", best_label, "NA"))
        elif second_score is not None and (best_score - second_score) < args.margin:
            rows.append((unit_id, best_label, "ambiguous_nearest", f"{best_score:.4f}", second_label, f"{second_score:.4f}"))
        else:
            rows.append((unit_id, best_label, "nearest_anchor", f"{best_score:.4f}", second_label or "NA", f"{second_score:.4f}" if second_score is not None else "NA"))

    with open(args.output, "w") as out_fh:
        out_fh.write("unit_id\tassigned_label\tmethod\tbest_avg_peq\tsecond_best_label\tsecond_best_avg_peq\n")
        for row in sorted(rows):
            out_fh.write("\t".join(row) + "\n")

    n_anchor = sum(1 for r in rows if r[2] == "anchor")
    n_nearest = sum(1 for r in rows if r[2] == "nearest_anchor")
    n_ambiguous = sum(1 for r in rows if r[2] == "ambiguous_nearest")
    n_unclassified = sum(1 for r in rows if r[2] in ("below_floor", "no_anchor_data"))
    n_conflict = sum(1 for r in rows if r[2] == "anchor_conflict")

    print(f"Total units: {len(rows)}")
    print(f"  anchor (ground truth): {n_anchor}")
    print(f"  nearest_anchor (assigned by similarity): {n_nearest}")
    print(f"  ambiguous_nearest (close call, flagged): {n_ambiguous}")
    print(f"  unclassified (below floor or no data): {n_unclassified}")
    if n_conflict:
        print(f"  anchor_conflict (unit matched >1 paper cluster in crosswalk - check manually): {n_conflict}")
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
