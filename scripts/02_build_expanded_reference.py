import argparse
import hashlib
from pathlib import Path
from Bio import SeqIO

def seq_hash(seq):
    return hashlib.sha256(str(seq).upper().encode()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--existing-reference", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-length", type=int, default=5000)
    ap.add_argument("--max-length", type=int, default=300000)
    args = ap.parse_args()

    plasmid_dir = Path(args.plasmid_fastas_dir)
    seen_hashes = {}
    records_out = []
    n_from_existing = 0
    n_from_cohort = 0
    n_duplicates_skipped = 0
    n_too_short_skipped = 0
    n_too_long_skipped = 0

    if args.existing_reference:
        for r in SeqIO.parse(args.existing_reference, "fasta"):
            if len(r.seq) > args.max_length:
                n_too_long_skipped += 1
                continue
            h = seq_hash(r.seq)
            if h in seen_hashes:
                n_duplicates_skipped += 1
                continue
            seen_hashes[h] = r.id
            records_out.append(r)
            n_from_existing += 1

    for fasta_path in sorted(plasmid_dir.glob("*_plasmids.fasta")):
        genome_id = fasta_path.stem.replace("_plasmids", "")
        for r in SeqIO.parse(fasta_path, "fasta"):
            if len(r.seq) < args.min_length:
                n_too_short_skipped += 1
                continue
            if len(r.seq) > args.max_length:
                n_too_long_skipped += 1
                continue
            h = seq_hash(r.seq)
            if h in seen_hashes:
                n_duplicates_skipped += 1
                continue
            seen_hashes[h] = r.id
            if not r.id.startswith(genome_id):
                r.id = f"{genome_id}__{r.id}"
            r.description = f"length_{len(r.seq)}bp source_cohort"
            records_out.append(r)
            n_from_cohort += 1

    with open(args.output, "w") as out_fh:
        SeqIO.write(records_out, out_fh, "fasta")

    print(f"Wrote {len(records_out)} unique plasmid sequences to {args.output}")
    print(f"  from existing reference: {n_from_existing}")
    print(f"  from cohort-resolved plasmids: {n_from_cohort}")
    print(f"  duplicate sequences skipped: {n_duplicates_skipped}")
    print(f"  too-short sequences skipped (<{args.min_length}bp): {n_too_short_skipped}")
    print(f"  too-long sequences skipped (>{args.max_length}bp): {n_too_long_skipped}")

if __name__ == "__main__":
    main()
