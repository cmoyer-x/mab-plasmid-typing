import argparse
import csv
from pathlib import Path
from Bio import SeqIO

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
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--all-genes-faa", required=True)
    ap.add_argument("--type-label", required=True, help="e.g. pB, pG, singleton_pGD13")
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()

    genome_ids = load_genome_ids(args.summary_tsv)

    member_units = []
    with open(args.plasmid_types_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["final_type_label"] == args.type_label:
                member_units.append(row["unit_id"])

    if not member_units:
        raise SystemExit(f"No units found with final_type_label == {args.type_label!r}")

    print(f"{len(member_units)} plasmid units belong to type {args.type_label!r}")

    plasmid_dir = Path(args.plasmid_fastas_dir)
    nuc_records = []
    for unit_id in member_units:
        gid = resolve_genome_id(unit_id, genome_ids)
        if gid is None:
            print(f"WARNING: could not resolve genome for {unit_id}, skipping")
            continue
        fasta_path = plasmid_dir / f"{gid}_plasmids.fasta"
        if not fasta_path.exists():
            print(f"WARNING: {fasta_path} not found, skipping {unit_id}")
            continue
        for r in SeqIO.parse(fasta_path, "fasta"):
            candidate_id = r.id if r.id.startswith(gid) else f"{gid}__{r.id}"
            if candidate_id == unit_id:
                nuc_records.append(r)
                break

    nuc_out = f"{args.output_prefix}_nucleotide.fasta"
    SeqIO.write(nuc_records, nuc_out, "fasta")
    print(f"Wrote {len(nuc_records)} nucleotide sequences to {nuc_out}")

    prot_records = []
    for r in SeqIO.parse(args.all_genes_faa, "fasta"):
        gene_unit_id = r.id.rsplit("|gene", 1)[0]
        if gene_unit_id in member_units:
            prot_records.append(r)

    prot_out = f"{args.output_prefix}_proteins.faa"
    SeqIO.write(prot_records, prot_out, "fasta")
    print(f"Wrote {len(prot_records)} protein sequences to {prot_out}")

if __name__ == "__main__":
    main()
