import argparse
import re
import subprocess
from pathlib import Path
from Bio import SeqIO

def find_circular(records, max_length):
    keep = []
    excluded = []
    for r in records:
        if "circular=true" in r.description.lower():
            if len(r.seq) <= max_length:
                keep.append(r)
            else:
                excluded.append(r)
    return keep, excluded

def find_handcurated(records, pattern):
    rx = re.compile(pattern, re.IGNORECASE)
    return [r for r in records if rx.search(r.id)]

def dedupe(records):
    seen = {}
    for r in records:
        seen[r.id] = r
    return list(seen.values())

def write_fasta(records, path):
    if records:
        SeqIO.write(records, path, "fasta")

def parse_reference_strain_map(reference_plasmids):
    strain_map = {}
    rx = re.compile(r"strain_(GD\d+[A-Za-z]*)", re.IGNORECASE)
    for r in SeqIO.parse(reference_plasmids, "fasta"):
        m = rx.search(r.description)
        if m:
            strain_map.setdefault(m.group(1).upper(), []).append(r)
    return strain_map

def run_blast(query_fasta, db_fasta, out_tsv, pident_min, cov_min, work_prefix):
    db_path = work_prefix.with_suffix("")
    subprocess.run(
        ["makeblastdb", "-in", str(db_fasta), "-dbtype", "nucl", "-out", str(db_path)],
        check=True
    )
    subprocess.run(
        [
            "blastn",
            "-query", str(query_fasta),
            "-db", str(db_path),
            "-out", str(out_tsv),
            "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore qcovs",
            "-evalue", "1e-10"
        ],
        check=True
    )
    best_qcovs = {}
    best_pident = {}
    if out_tsv.exists():
        with open(out_tsv) as fh:
            for line in fh:
                fields = line.strip().split("\t")
                qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore, qcovs = fields
                pident = float(pident)
                qcovs = float(qcovs) / 100.0
                key = (qseqid, sseqid)
                if qcovs > best_qcovs.get(key, 0):
                    best_qcovs[key] = qcovs
                    best_pident[key] = pident

    hits = {}
    for (qseqid, sseqid), qcovs in best_qcovs.items():
        pident = best_pident[(qseqid, sseqid)]
        if pident >= pident_min and qcovs >= cov_min:
            hits.setdefault(qseqid, []).append((sseqid, pident, qcovs))
    return hits

def screen_phix(records, phix_reference, blast_dir, pident_min, cov_min):
    if not records or phix_reference is None:
        return records, []

    query_fasta = blast_dir / "__phix_screen_query.fasta"
    write_fasta(records, query_fasta)

    out_tsv = blast_dir / "__phix_screen_results.tsv"
    hits = run_blast(
        query_fasta, phix_reference, out_tsv,
        pident_min, cov_min, blast_dir / "__phix_screen_db"
    )

    clean = [r for r in records if r.id not in hits]
    contaminated = [r for r in records if r.id in hits]
    return clean, contaminated

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemblies-dir", required=True)
    ap.add_argument("--reference-plasmids", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--handcurated-pattern", default=r"pGD\d+")
    ap.add_argument("--pident-min", type=float, default=75.0)
    ap.add_argument("--cov-min", type=float, default=0.10)
    ap.add_argument("--direct-pident-min", type=float, default=98.0)
    ap.add_argument("--direct-cov-min", type=float, default=0.9)
    ap.add_argument("--max-plasmid-length", type=int, default=300000)
    ap.add_argument("--phix-reference", default=None)
    ap.add_argument("--phix-pident-min", type=float, default=90.0)
    ap.add_argument("--phix-cov-min", type=float, default=0.5)
    args = ap.parse_args()

    assemblies_dir = Path(args.assemblies_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plasmid_dir = output_dir / "plasmid_fastas"
    plasmid_dir.mkdir(exist_ok=True)
    blast_dir = output_dir / "blast_tmp"
    blast_dir.mkdir(exist_ok=True)

    phix_reference = Path(args.phix_reference) if args.phix_reference else None
    phix_contaminants = []

    reference_plasmids = Path(args.reference_plasmids)
    strain_map = parse_reference_strain_map(reference_plasmids)

    summary_rows = []
    excluded_rows = []
    pending_direct = {}
    pending_pooled = {}

    for fasta_path in sorted(assemblies_dir.glob("*.fasta")):
        genome_id = fasta_path.stem
        records = list(SeqIO.parse(fasta_path, "fasta"))

        circular_hits, excluded_circular = find_circular(records, args.max_plasmid_length)
        handcurated_hits = find_handcurated(records, args.handcurated_pattern)
        combined = dedupe(circular_hits + handcurated_hits)

        combined, phix_hits = screen_phix(
            combined, phix_reference, blast_dir, args.phix_pident_min, args.phix_cov_min
        )
        for r in phix_hits:
            phix_contaminants.append((genome_id, r.id, len(r.seq)))

        for r in excluded_circular:
            excluded_rows.append((genome_id, r.id, len(r.seq)))

        if combined:
            out_path = plasmid_dir / f"{genome_id}_plasmids.fasta"
            write_fasta(combined, out_path)
            method = []
            if circular_hits:
                method.append("circular")
            if handcurated_hits:
                method.append("handcurated")
            summary_rows.append([genome_id, "+".join(method), ",".join(r.id for r in combined)])
            continue

        strain_key = None
        for key in strain_map:
            if genome_id.upper().startswith(key):
                strain_key = key
                break

        if strain_key:
            pending_direct[genome_id] = (records, strain_map[strain_key])
        else:
            pending_pooled[genome_id] = records

        summary_rows.append([genome_id, "pending", ""])

    if pending_direct:
        direct_query = blast_dir / "direct_match_contigs.fasta"
        direct_ref = blast_dir / "direct_match_reference.fasta"
        with open(direct_query, "w") as qfh, open(direct_ref, "w") as rfh:
            written_refs = set()
            for genome_id, (records, ref_records) in pending_direct.items():
                for r in records:
                    r.id = f"{genome_id}__{r.id}"
                    r.description = ""
                SeqIO.write(records, qfh, "fasta")
                for ref in ref_records:
                    if ref.id not in written_refs:
                        SeqIO.write([ref], rfh, "fasta")
                        written_refs.add(ref.id)

        direct_out = blast_dir / "direct_match_results.tsv"
        direct_hits = run_blast(
            direct_query, direct_ref, direct_out,
            args.direct_pident_min, args.direct_cov_min,
            blast_dir / "direct_match_db"
        )

        for i, row in enumerate(summary_rows):
            genome_id = row[0]
            if genome_id not in pending_direct or row[1] != "pending":
                continue
            records, ref_records = pending_direct[genome_id]
            matched_contigs = []
            for r in records:
                if r.id in direct_hits:
                    matched_contigs.append(r)
            if matched_contigs:
                out_path = plasmid_dir / f"{genome_id}_plasmids.fasta"
                write_fasta(matched_contigs, out_path)
                summary_rows[i] = [genome_id, "direct_strain_reference", ",".join(r.id for r in matched_contigs)]
            else:
                pending_pooled[genome_id] = records
                summary_rows[i] = [genome_id, "pending", ""]

    if pending_pooled:
        combined_query = blast_dir / "pooled_contigs.fasta"
        with open(combined_query, "w") as out_fh:
            for genome_id, records in pending_pooled.items():
                for r in records:
                    r.id = f"{genome_id}__{r.id}"
                    r.description = ""
                SeqIO.write(records, out_fh, "fasta")

        pooled_out = blast_dir / "pooled_results.tsv"
        pooled_hits = run_blast(
            combined_query, reference_plasmids, pooled_out,
            args.pident_min, args.cov_min,
            blast_dir / "pooled_db"
        )

        hits_by_genome = {}
        for hit_id in pooled_hits:
            genome_id, contig_id = hit_id.split("__", 1)
            hits_by_genome.setdefault(genome_id, []).append(hit_id)

        for i, row in enumerate(summary_rows):
            genome_id = row[0]
            if genome_id not in pending_pooled or row[1] != "pending":
                continue
            records = pending_pooled[genome_id]
            matched_ids = set(hits_by_genome.get(genome_id, []))
            matched_records = [
                r for r in records
                if r.id in matched_ids and len(r.seq) <= args.max_plasmid_length
            ]
            if matched_records:
                out_path = plasmid_dir / f"{genome_id}_plasmids.fasta"
                write_fasta(matched_records, out_path)
                summary_rows[i] = [genome_id, "blast_homology", ",".join(r.id for r in matched_records)]
            else:
                summary_rows[i] = [genome_id, "no_plasmid_detected", ""]

    summary_path = output_dir / "plasmid_extraction_summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write("genome_id\tmethod\tcontig_ids\n")
        for genome_id, method, ids in summary_rows:
            fh.write(f"{genome_id}\t{method}\t{ids}\n")

    excluded_path = output_dir / "excluded_full_length_circular_contigs.tsv"
    with open(excluded_path, "w") as fh:
        fh.write("genome_id\tcontig_id\tlength\n")
        for genome_id, contig_id, length in excluded_rows:
            fh.write(f"{genome_id}\t{contig_id}\t{length}\n")

    phix_path = output_dir / "excluded_phix_contamination.tsv"
    with open(phix_path, "w") as fh:
        fh.write("genome_id\tcontig_id\tlength\n")
        for genome_id, contig_id, length in phix_contaminants:
            fh.write(f"{genome_id}\t{contig_id}\t{length}\n")

    print(f"Summary written to {summary_path}")
    print(f"Excluded full-length circular contigs written to {excluded_path}")
    if phix_reference is not None:
        print(f"Excluded PhiX174 contamination ({len(phix_contaminants)} contigs) written to {phix_path}")
    print(f"Plasmid FASTAs written to {plasmid_dir}")

if __name__ == "__main__":
    main()
