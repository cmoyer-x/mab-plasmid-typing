import argparse
from pathlib import Path
from Bio import SeqIO

def load_gene_to_pham(pham_raw_dir):
    gene_to_pham = {}
    for pham_fasta in sorted(Path(pham_raw_dir).rglob("*.faa")):
        pham_id = pham_fasta.stem
        for r in SeqIO.parse(pham_fasta, "fasta"):
            gene_to_pham[r.id] = pham_id
    return gene_to_pham

def unit_phams(unit_id, all_genes_faa, gene_to_pham):
    phams = set()
    for r in SeqIO.parse(all_genes_faa, "fasta"):
        if r.id.startswith(f"{unit_id}|gene"):
            pham_id = gene_to_pham.get(r.id, f"UNASSIGNED_{r.id}")
            phams.add(pham_id)
    return phams

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pham-raw-dir", required=True)
    ap.add_argument("--all-genes-faa", required=True)
    ap.add_argument("--pairs", required=True, nargs="+",
                     help="pairs as smaller:larger,smaller:larger,...")
    args = ap.parse_args()

    gene_to_pham = load_gene_to_pham(args.pham_raw_dir)

    for pair_str in args.pairs:
        smaller, larger = pair_str.split(":")
        smaller_phams = unit_phams(smaller, args.all_genes_faa, gene_to_pham)
        larger_phams = unit_phams(larger, args.all_genes_faa, gene_to_pham)

        shared = smaller_phams & larger_phams
        only_smaller = smaller_phams - larger_phams
        only_larger = larger_phams - smaller_phams

        pct_smaller_contained = 100 * len(shared) / len(smaller_phams) if smaller_phams else 0

        print(f"=== {smaller} ({len(smaller_phams)} phams) vs {larger} ({len(larger_phams)} phams) ===")
        print(f"  shared phams: {len(shared)}")
        print(f"  phams only in {smaller}: {len(only_smaller)}")
        print(f"  phams only in {larger}: {len(only_larger)}")
        print(f"  % of {smaller}'s genes contained in {larger}: {pct_smaller_contained:.1f}%")
        if only_smaller:
            print(f"  {smaller}-only phams: {sorted(only_smaller)}")
        print()

if __name__ == "__main__":
    main()
