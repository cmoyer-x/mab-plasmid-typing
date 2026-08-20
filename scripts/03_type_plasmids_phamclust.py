import argparse
import csv
import subprocess
from pathlib import Path
from Bio import SeqIO

COMPLETE_METHODS = {"circular", "handcurated", "direct_strain_reference", "circular+handcurated"}

def load_summary(summary_tsv):
    rows = {}
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows[row["genome_id"]] = row
    return rows

def gather_complete_units(plasmid_fastas_dir, summary):
    units = {}
    for genome_id, row in summary.items():
        if row["method"] not in COMPLETE_METHODS:
            continue
        fasta_path = plasmid_fastas_dir / f"{genome_id}_plasmids.fasta"
        if not fasta_path.exists():
            continue
        for r in SeqIO.parse(fasta_path, "fasta"):
            unit_id = r.id if r.id.startswith(genome_id) else f"{genome_id}__{r.id}"
            units[unit_id] = r.seq
    return units

def run_prodigal(unit_id, seq, work_dir):
    nuc_fasta = work_dir / f"{unit_id}.fna"
    with open(nuc_fasta, "w") as fh:
        fh.write(f">{unit_id}\n{str(seq)}\n")

    faa_out = work_dir / f"{unit_id}.faa"
    subprocess.run(
        ["prodigal", "-i", str(nuc_fasta), "-a", str(faa_out), "-p", "meta", "-q"],
        check=True
    )

    genes = []
    for i, r in enumerate(SeqIO.parse(faa_out, "fasta"), start=1):
        genes.append((f"{unit_id}|gene{i}", str(r.seq).rstrip("*")))
    return genes

def run(cmd):
    subprocess.run(cmd, check=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--clu-thresh", type=float, default=None)
    args = ap.parse_args()

    plasmid_dir = Path(args.plasmid_fastas_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genecall_dir = output_dir / "genecalling"
    genecall_dir.mkdir(exist_ok=True)

    summary = load_summary(args.summary_tsv)
    units = gather_complete_units(plasmid_dir, summary)
    print(f"Found {len(units)} complete plasmid units to type")

    all_genes_faa = output_dir / "all_genes.faa"
    unit_gene_lists = {}

    if all_genes_faa.exists() and not args.force:
        print(f"{all_genes_faa} already exists, skipping gene calling (use --force to redo)")
        current_unit_id = None
        for r in SeqIO.parse(all_genes_faa, "fasta"):
            unit_id = r.id.rsplit("|gene", 1)[0]
            unit_gene_lists.setdefault(unit_id, []).append(r.id)
    else:
        with open(all_genes_faa, "w") as out_fh:
            for i, (unit_id, seq) in enumerate(units.items(), start=1):
                genes = run_prodigal(unit_id, seq, genecall_dir)
                unit_gene_lists[unit_id] = [g[0] for g in genes]
                for gene_id, translation in genes:
                    out_fh.write(f">{gene_id}\n{translation}\n")
                if i % 25 == 0:
                    print(f"  called genes for {i}/{len(units)} units")

    total_genes = sum(len(g) for g in unit_gene_lists.values())
    print(f"Have {total_genes} genes across {len(unit_gene_lists)} plasmid units")

    pham_raw_dir = output_dir / "pham_raw"
    if pham_raw_dir.exists() and not args.force:
        print(f"{pham_raw_dir} already exists, skipping phammseqs (use --force to redo)")
    else:
        run([
            "phammseqs", str(all_genes_faa),
            "-o", str(pham_raw_dir),
            "-c", str(args.threads),
            "-v"
        ])

    gene_to_pham = {}
    for pham_fasta in sorted(pham_raw_dir.rglob("*.faa")) + sorted(pham_raw_dir.rglob("*.fasta")) + sorted(pham_raw_dir.rglob("*.fa")):
        pham_id = pham_fasta.stem
        for r in SeqIO.parse(pham_fasta, "fasta"):
            gene_to_pham[r.id] = pham_id

    unassigned = [g for genes in unit_gene_lists.values() for g in genes if g not in gene_to_pham]
    if unassigned:
        print(f"WARNING: {len(unassigned)} genes were not found in any pham output file")
        print(f"  example: {unassigned[0]}")

    phamclust_input_dir = output_dir / "phamclust_input"
    phamclust_input_dir.mkdir(exist_ok=True)

    all_genes_by_id = {}
    for r in SeqIO.parse(all_genes_faa, "fasta"):
        all_genes_by_id[r.id] = str(r.seq)

    for unit_id, gene_ids in unit_gene_lists.items():
        pham_counts = {}
        out_path = phamclust_input_dir / f"{unit_id}.faa"
        with open(out_path, "w") as fh:
            for gene_id in gene_ids:
                pham_id = gene_to_pham.get(gene_id)
                if pham_id is None:
                    continue
                pham_counts[pham_id] = pham_counts.get(pham_id, 0) + 1
                n = pham_counts[pham_id]
                translation = all_genes_by_id[gene_id]
                fh.write(f">name={unit_id}|pham={pham_id}|n={n}\n{translation}\n")

    phamclust_out_dir = output_dir / "phamclust_results"
    phamclust_cmd = [
        "phamclust", str(phamclust_input_dir), str(phamclust_out_dir),
        "-g", "-t", str(args.threads)
    ]
    if args.clu_thresh is not None:
        phamclust_cmd += ["-c", str(args.clu_thresh)]
    run(phamclust_cmd)

    print(f"phamclust results written to {phamclust_out_dir}")

if __name__ == "__main__":
    main()
