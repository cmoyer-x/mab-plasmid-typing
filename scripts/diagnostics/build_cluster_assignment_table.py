import argparse
import re
from pathlib import Path

def read_matrix_ids(tsv_path):
    with open(tsv_path) as fh:
        first_line = fh.readline().rstrip("\n")
    fields = first_line.split("\t")
    ids = [f for f in fields if f]
    return ids

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phamclust-results-dir", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    results_dir = Path(args.phamclust_results_dir)
    assignments = {}

    cluster_dirs = sorted(
        [d for d in results_dir.iterdir() if d.is_dir() and re.match(r"cluster_\d+$", d.name)],
        key=lambda d: int(d.name.split("_")[1])
    )

    for cluster_dir in cluster_dirs:
        cluster_id = cluster_dir.name
        peq_file = cluster_dir / "peq_similarity.tsv"
        if peq_file.exists():
            for unit_id in read_matrix_ids(peq_file):
                assignments[unit_id] = {"cluster": cluster_id, "subcluster": "NA"}

        subcluster_files = sorted(
            cluster_dir.glob("subcluster_*_similarity.tsv"),
            key=lambda p: int(re.search(r"subcluster_(\d+)_", p.name).group(1))
        )
        for subcluster_file in subcluster_files:
            subcluster_num = re.search(r"subcluster_(\d+)_", subcluster_file.name).group(1)
            subcluster_id = f"{cluster_id}_sub{subcluster_num}"
            for unit_id in read_matrix_ids(subcluster_file):
                if unit_id in assignments and assignments[unit_id]["cluster"] == cluster_id:
                    assignments[unit_id]["subcluster"] = subcluster_id

    singletons_dir = results_dir / "singletons" / "genomes"
    if singletons_dir.exists():
        for f in singletons_dir.glob("*.faa"):
            unit_id = f.stem
            assignments[unit_id] = {"cluster": "singleton", "subcluster": "NA"}

    with open(args.output, "w") as out_fh:
        out_fh.write("unit_id\tcluster_id\tsubcluster_id\n")
        for unit_id, info in sorted(assignments.items()):
            out_fh.write(f"{unit_id}\t{info['cluster']}\t{info['subcluster']}\n")

    n_clusters = len(cluster_dirs)
    n_singletons = sum(1 for v in assignments.values() if v["cluster"] == "singleton")
    print(f"Wrote {len(assignments)} genome-to-cluster assignments to {args.output}")
    print(f"  {n_clusters} clusters, {n_singletons} singletons")

if __name__ == "__main__":
    main()
