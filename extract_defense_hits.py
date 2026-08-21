import argparse
import csv

def parse_fasta(path):
    records = {}
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records[header] = "".join(seq_lines)
    return records

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes-tsv", required=True)
    ap.add_argument("--prt-fasta", required=True)
    ap.add_argument("--gene-name", default="VP1853__VP1853")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    hit_ids = []
    with open(args.genes_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["gene_name"] == args.gene_name:
                hit_ids.append((row["replicon"], row["hit_id"]))

    print(f"{len(hit_ids)} hits found for {args.gene_name} in {args.genes_tsv}")

    proteins = parse_fasta(args.prt_fasta)

    written = 0
    missing = []
    with open(args.output, "w") as out_fh:
        for replicon, hit_id in hit_ids:
            if hit_id not in proteins:
                missing.append(hit_id)
                continue
            seq = proteins[hit_id]
            out_fh.write(f">{hit_id} replicon={replicon} length={len(seq)}aa\n")
            for i in range(0, len(seq), 70):
                out_fh.write(seq[i:i+70] + "\n")
            written += 1

    print(f"Wrote {written} protein sequences to {args.output}")
    if missing:
        print(f"WARNING: {len(missing)} hit_ids not found in {args.prt_fasta}: {missing}")

if __name__ == "__main__":
    main()
