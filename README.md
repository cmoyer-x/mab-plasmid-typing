# M. abscessus cohort plasmid extraction and typing

Extracts, expands, and types plasmids across a 425-genome *M. abscessus*
clinical cohort, cross-referenced against the plasmid clusters defined in
Dedrick et al. (mBio, 2021), *"The Prophage and Plasmid Mobilome as a Likely
Driver of Mycobacterium abscessus Diversity."*

## Pipeline overview

```
01_extract_plasmids.py            per-genome plasmid extraction (circular/
                                   hand-curated/direct-match/BLAST homology),
                                   with PhiX174 contamination screening
02_build_expanded_reference.py    pools resolved plasmids into an expanded
                                   BLAST reference to catch within-cohort HGT
03_type_plasmids_phamclust.py     gene-content typing (Prodigal -> phammseqs
                                   -> phamclust) producing a PEQ similarity
                                   matrix
04_assign_final_plasmid_types.py  final typing: matches genomes to the
                                   paper's clusters via exact-length anchors,
                                   assigns everything else by nearest-anchor
                                   similarity, and clusters genuinely novel
                                   plasmids into new pI/pJ/pK... types
```

`scripts/utils/compare_gene_content.py` checks whether one plasmid's gene
content is a subset of another's — used to confirm a backbone+cargo
relationship between two of the novel clusters.

`scripts/diagnostics/` holds the exploratory scripts used to validate the
methodology (phamclust cluster crosswalk, nearest-anchor assignment,
unclassified-genome characterization). Their logic was folded into
`04_assign_final_plasmid_types.py`; kept here for transparency into how that
script's approach was derived and validated.

## Dependencies

```
conda install -c bioconda blast prodigal mmseqs2=13.45111 clustalo -y
pip install biopython phammseqs phamclust
```

## Usage

```bash
# 1. Extract plasmids per genome
python scripts/01_extract_plasmids.py \
  --assemblies-dir /path/to/ALL_GD_fasta \
  --reference-plasmids reference/reference_plasmids_Mabscessus.fasta \
  --phix-reference reference/phix174_reference.fasta \
  --output-dir plasmids/

# 2. Expand the reference with cohort-resolved plasmids (repeat as useful)
python scripts/02_build_expanded_reference.py \
  --plasmid-fastas-dir plasmids/plasmid_fastas \
  --existing-reference reference/reference_plasmids_Mabscessus.fasta \
  --output reference_plasmids_expanded.fasta

# 3. Re-run extraction against the expanded reference
python scripts/01_extract_plasmids.py \
  --assemblies-dir /path/to/ALL_GD_fasta \
  --reference-plasmids reference_plasmids_expanded.fasta \
  --phix-reference reference/phix174_reference.fasta \
  --output-dir plasmids_expanded/

# 4. Gene-content typing
python scripts/03_type_plasmids_phamclust.py \
  --plasmid-fastas-dir plasmids_expanded/plasmid_fastas \
  --summary-tsv plasmids_expanded/plasmid_extraction_summary.tsv \
  --output-dir plasmid_phamclust/ \
  --threads 8

# 5. Final typing (paper clusters + novel pI/pJ/pK... clusters)
python scripts/04_assign_final_plasmid_types.py \
  --summary-tsv plasmids_expanded/plasmid_extraction_summary.tsv \
  --plasmid-fastas-dir plasmids_expanded/plasmid_fastas \
  --all-genes-faa plasmid_phamclust/all_genes.faa \
  --similarity-matrix plasmid_phamclust/phamclust_results/pairwise_peq_similarities.tsv \
  --output plasmid_final_types.tsv
```

## Methodology notes / lessons learned

- **Whole-sequence Mash distance clustering was tried first and abandoned.**
  It couldn't distinguish plasmids that share one conserved replication gene
  (e.g. the paper's proposed "IncMabI" group: pA/pC/pD/pE/pF) from plasmids
  that are actually the same type. Gene-content clustering (phammseqs/
  phamclust) resolves this because it isn't fooled by a single shared gene.
- **phamclust's own cluster boundaries have linkage-sensitivity artifacts**
  in a heterogeneous set — verified by checking specific genomes (e.g. one
  unit was split into its own cluster despite scoring *higher* similarity to
  the cluster it left than to anything else). `04_assign_final_plasmid_types.py`
  works around this by using paper-verified anchors plus direct nearest-
  neighbor similarity rather than trusting phamclust's cluster assignment
  directly.
- **PhiX174 contamination**: three genomes' "novel" plasmid cluster turned
  out to be the Illumina sequencing control genome (exactly 5,386 bp,
  100% match to NC_001422.1) that wasn't filtered from reads before
  assembly. Both `01_extract_plasmids.py` (via `--phix-reference`) and
  `04_assign_final_plasmid_types.py` (via exact-length check) now screen
  for this.
- **Two genomes matching `pATCC19977`'s exact length (23,319 bp) currently
  sit in the wrong novel cluster** because the external `pMAB23` reference
  sequence was never run through gene-calling, so no anchor exists for it.
  Fix: add `pMAB23`'s sequence into the `03_type_plasmids_phamclust.py`
  gene-calling input.

## Current cohort result

- 425 genomes total; 347 with an identified plasmid before typing (81.6%).
- 210 complete plasmid units typed (after excluding PhiX contamination and
  short fragments): 154 match a paper-defined cluster (pA-pH or a named
  singleton), 51 form 10 novel cohort-specific clusters (pI-pR) plus
  singletons.
- Notable novel finding: a 10-member ~13.5 kb plasmid type (pI) with no
  representative in the original 2021 catalog.
- Notable relationship: two novel clusters each have a larger relative whose
  gene content is a near-total superset (95-100%) of the smaller type plus
  16-27 additional gene families - consistent with a fusion/cointegrate or
  deletion-derivative relationship, worth dedicated follow-up.
