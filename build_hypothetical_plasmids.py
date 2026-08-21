import argparse
import csv
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def load_best_reference(pooled_results_tsv):
    best_score = {}
    best_ref = {}
    with open(pooled_results_tsv) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore, qcovs = fields
            qcovs = float(qcovs)
            pident = float(pident)
            score = (qcovs, pident)
            if qseqid not in best_score or score > best_score[qseqid]:
                best_score[qseqid] = score
                best_ref[qseqid] = sseqid
    return best_ref

def load_genome_ids(summary_tsv):
    genome_ids = set()
    blast_homology_contigs = {}
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_ids.add(row["genome_id"])
            if row["method"] == "blast_homology" and row["contig_ids"]:
                blast_homology_contigs[row["genome_id"]] = row["contig_ids"].split(",")
    return genome_ids, blast_homology_contigs

def load_plasmid_records(plasmid_fastas_dir, genome_id):
    fasta_path = Path(plasmid_fastas_dir) / f"{genome_id}_plasmids.fasta"
    records = {}
    if fasta_path.exists():
        for r in SeqIO.parse(fasta_path, "fasta"):
            records[r.id] = r
    return records

def load_reference_lengths(reference_fasta):
    lengths = {}
    header = None
    seq_len = 0
    with open(reference_fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    lengths[header] = seq_len
                header = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
        if header is not None:
            lengths[header] = seq_len
    return lengths

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--pooled-results-tsv", required=True)
    ap.add_argument("--plasmid-fastas-dir", required=True)
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--spacer-length", type=int, default=100)
    ap.add_argument("--max-length-ratio", type=float, default=1.5,
                     help="max allowed ratio of combined length to reference length")
    ap.add_argument("--max-plasmid-length", type=int, default=300000)
    args = ap.parse_args()

    genome_ids, blast_homology_contigs = load_genome_ids(args.summary_tsv)
    best_ref = load_best_reference(args.pooled_results_tsv)
    reference_lengths = load_reference_lengths(args.reference_fasta)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    hypothetical_dir = output_dir / "hypothetical_plasmid_fastas"
    hypothetical_dir.mkdir(exist_ok=True)

    manifest_rows = []
    rejected_rows = []
    spacer = "N" * args.spacer_length

    for genome_id, contig_ids in blast_homology_contigs.items():
        records = load_plasmid_records(args.plasmid_fastas_dir, genome_id)

        groups = defaultdict(list)
        for contig_id in contig_ids:
            qseqid = f"{genome_id}__{contig_id}" if not contig_id.startswith(genome_id) else contig_id
            ref = best_ref.get(qseqid, "UNKNOWN_REFERENCE")
            groups[ref].append(contig_id)

        out_records = []
        hyp_counter = 1
        for ref, member_contigs in groups.items():
            if len(member_contigs) < 2:
                continue

            seqs = []
            for contig_id in member_contigs:
                rec = records.get(contig_id)
                if rec is None:
                    continue
                seqs.append(str(rec.seq))

            if len(seqs) < 2:
                continue

            combined_seq = spacer.join(seqs)
            combined_length = len(combined_seq)
            ref_length = reference_lengths.get(ref)

            if combined_length > args.max_plasmid_length:
                rejected_rows.append((
                    genome_id, ref, len(member_contigs), ",".join(member_contigs),
                    combined_length, ref_length, "exceeds_max_plasmid_length"
                ))
                continue

            if ref_length is not None:
                ratio = combined_length / ref_length
                if ratio > args.max_length_ratio:
                    rejected_rows.append((
                        genome_id, ref, len(member_contigs), ",".join(member_contigs),
                        combined_length, ref_length, f"length_ratio_{ratio:.1f}x_exceeds_{args.max_length_ratio}x"
                    ))
                    continue
            else:
                rejected_rows.append((
                    genome_id, ref, len(member_contigs), ",".join(member_contigs),
                    combined_length, "NA", "reference_length_unknown"
                ))
                continue

            hyp_id = f"{genome_id}_hypothetical{hyp_counter}"
            hyp_counter += 1

            new_record = SeqRecord(
                Seq(combined_seq),
                id=hyp_id,
                description=(
                    f"hypothetical_plasmid n_fragments={len(seqs)} "
                    f"reference={ref} length={combined_length}bp "
                    f"members={','.join(member_contigs)}"
                )
            )
            out_records.append(new_record)

            manifest_rows.append((
                hyp_id, genome_id, ref, len(member_contigs),
                ",".join(member_contigs), combined_length
            ))

        if out_records:
            out_path = hypothetical_dir / f"{genome_id}_hypothetical.fasta"
            SeqIO.write(out_records, out_path, "fasta")

    manifest_path = output_dir / "hypothetical_plasmids_manifest.tsv"
    with open(manifest_path, "w") as fh:
        fh.write("hypothetical_id\tgenome_id\treference_subject\tn_fragments\tmember_contig_ids\ttotal_length_bp\n")
        for row in manifest_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    rejected_path = output_dir / "hypothetical_plasmids_rejected.tsv"
    with open(rejected_path, "w") as fh:
        fh.write("genome_id\treference_subject\tn_fragments\tmember_contig_ids\tcombined_length_bp\treference_length_bp\treject_reason\n")
        for row in rejected_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{len(manifest_rows)} hypothetical plasmids built and accepted")
    print(f"{len(rejected_rows)} groups rejected as implausible (see {rejected_path})")
    print(f"Manifest written to {manifest_path}")
    print(f"FASTAs written to {hypothetical_dir}")

if __name__ == "__main__":
    main()
