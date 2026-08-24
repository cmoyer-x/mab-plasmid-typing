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
05_build_hypothetical_plasmids.py combines multi-fragment BLAST-homology
                                   hits into "hypothetical" reconstructed
                                   plasmids when fragments share a common
                                   best-hit reference, with plausibility
                                   filtering against that reference's actual
                                   length to reject implausible/chimeric
                                   groupings (see methodology notes)
```

`scripts/analysis/` holds downstream analyses run on the final typing table:

- `correlate_plasmid_dcc.py` / `correlate_plasmid_eop.py` - test association
  between plasmid type and clonal complex (DCC) or phage susceptibility
  (EOP), with Fisher's exact / Mann-Whitney U and FDR correction
- `check_multi_plasmid_genomes.py` - identifies genomes carrying more than
  one distinct plasmid type
- `extract_type_sequences.py` - pulls nucleotide/protein sequences for a
  single plasmid type, ready for functional annotation (DefenseFinder,
  PADLOC, Bakta, etc.)
- `extract_defense_hits.py` - extracts specific DefenseFinder hit sequences
  from its raw `.prt` output by hit_id, for downstream alignment/comparison
- `extract_GD233_series.py` - pulls and compares plasmids across a specific
  longitudinal isolate series
- `plot_plasmid_dcc_heatmap.R` / `plot_cluster_sizes.R` - publication-style
  (pink/purple themed) visualizations of the association and clustering
  results

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

# 6. (optional) Reconstruct hypothetical plasmids from multi-fragment
#    BLAST-homology hits left over from step 3 (independent of steps 4-5;
#    only needs the extraction output, not the typing results)
python scripts/05_build_hypothetical_plasmids.py \
  --summary-tsv plasmids_expanded/plasmid_extraction_summary.tsv \
  --pooled-results-tsv plasmids_expanded/blast_tmp/pooled_results.tsv \
  --plasmid-fastas-dir plasmids_expanded/plasmid_fastas \
  --reference-fasta reference_plasmids_expanded.fasta \
  --output-dir hypothetical_plasmids/

# 7. (optional) Re-run gene-content typing including the hypothetical
#    plasmids - reuses already-called genes for the original complete
#    units, only calls genes for the new hypothetical ones, but always
#    rebuilds the pham/cluster space since new sequences change gene-family
#    assignments for everyone
python scripts/03_type_plasmids_phamclust.py \
  --plasmid-fastas-dir plasmids_expanded/plasmid_fastas \
  --summary-tsv plasmids_expanded/plasmid_extraction_summary.tsv \
  --hypothetical-fastas-dir hypothetical_plasmids/hypothetical_plasmid_fastas \
  --output-dir plasmid_phamclust_with_hypothetical/ \
  --threads 8

# 8. Final typing including hypothetical plasmids, with stable cluster
#    labels preserved from the previous (complete-units-only) run
python scripts/04_assign_final_plasmid_types.py \
  --summary-tsv plasmids_expanded/plasmid_extraction_summary.tsv \
  --plasmid-fastas-dir plasmids_expanded/plasmid_fastas \
  --hypothetical-fastas-dir hypothetical_plasmids/hypothetical_plasmid_fastas \
  --all-genes-faa plasmid_phamclust_with_hypothetical/all_genes.faa \
  --similarity-matrix plasmid_phamclust_with_hypothetical/phamclust_results/pairwise_peq_similarities.tsv \
  --previous-types-tsv plasmid_final_types.tsv \
  --output plasmid_final_types_with_hypothetical.tsv
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
- **Multi-fragment BLAST-homology hits require a plausibility check before
  combining into one hypothetical plasmid.** Grouping fragments purely by
  "same best-hit reference" is vulnerable to the same failure mode as the
  Mash/phamclust clustering issues above: a shared mobile element or
  conserved backbone gene can link genuinely unrelated fragments. Fix
  (`05_build_hypothetical_plasmids.py`): reject any combined group whose
  total length exceeds 300 kb outright, or exceeds 1.5x the actual
  reference plasmid's own length - this caught several 300-365 kb chimeric
  constructs built from an 80.8 kb reference (`pGD509`) before they could
  contaminate downstream typing.
- **Novel cluster letters (pI, pJ, pK...) are not stable across separate
  runs by default** - they get re-derived from scratch each run, ranked by
  cluster size in that specific run. Adding the 89 hypothetical plasmids
  changed the size ranking of every novel cluster and reshuffled most of
  the letters. Fix: `04_assign_final_plasmid_types.py` now accepts
  `--previous-types-tsv`, matching each new cluster against a prior run's
  labels by majority membership overlap (>=50%) and reusing that label
  when found, only minting a fresh letter for clusters with no real
  precedent. Always pass this pointing at the previous run's output when
  the underlying data changes (new genomes, added hypothetical plasmids,
  etc.) to keep type labels meaningful across the project's lifetime.

## Notable downstream findings

- **`pB` (35 carriers) is significantly associated with reduced phage
  susceptibility** (mean EOP fraction 0.0086 vs 0.0736 in non-carriers,
  FDR p = 0.00015) and **100% of pB plasmids carry a VP1853 antiphage
  defense system** (Doron et al. 2018 Science). `pB` is also enriched in
  genetically diverse `Non-DCC` backgrounds rather than one clone,
  consistent with repeated independent acquisition due to a fitness
  benefit.
- **`pG` (11 carriers, 100% restricted to DCC2) shows the opposite
  pattern**: significantly higher phage susceptibility (mean EOP 0.103 vs
  0.056, FDR p = 0.005) and 0% carry any detected defense system - a clean,
  mechanistically consistent counter-example.
- **`pL`, a novel cluster not in the original 2021 catalog, is perfectly
  restricted to DCC7** (4/4 carriers, FDR p = 5.3e-4) - strong evidence of
  clonal inheritance for a previously-undescribed plasmid type. (Note:
  earlier notes referred to this cluster as `pK` before the stable-labeling
  fix was in place - `pL` is the correct, stable label.)
- **`pC` (47 carriers, the largest type) is significantly excluded from
  DCC5** (0/47 carriers in that clonal complex) - a genuine depletion
  signal worth its own follow-up.


## Current cohort result

- 425 genomes total; 347 with an identified plasmid via extraction (81.6%).
- 317 plasmid units typed in total: 220 complete units (from `circular`,
  `handcurated`, `direct_strain_reference` methods) plus 97 hypothetical
  units reconstructed from multi-fragment BLAST-homology hits (89 of which
  passed the plausibility filter in `05_build_hypothetical_plasmids.py`;
  the rest were filtered as too short/too few genes once run through the
  same typing pipeline as everything else).
- Paper-cluster matches (anchor + nearest-anchor) and novel cohort-specific
  clusters (pI-pR-range letters, stably tracked across runs via
  `--previous-types-tsv`) span the full 317-unit set.
- Notable novel finding: a 10-member ~13.5 kb plasmid type (`pI`) with no
  representative in the original 2021 catalog.
- Notable relationship: `pK` (`GD10_13`/`GD272`/`GD276A-2`/`GD276B-2`, plus
  2 hypothetical members added after typing) has a smaller relative whose
  gene content is a near-total subset (95-100%) of the larger type - see
  methodology notes above (fusion/cointegrate or deletion-derivative
  relationship, worth dedicated follow-up).

## Open items

- [x] Type the hypothetical plasmids (done - see combined 317-unit result
  above).
- [x] Fix novel-cluster label stability across runs (done - see
  `--previous-types-tsv`).
- [ ] Add `pMAB23` into the gene-calling anchor set to resolve the
  `pATCC19977` mismatch noted above.
- [ ] Extend the DefenseFinder/EOP/DCC association checks to the remaining
  typed clusters (only pA, pB, pC, pD, pF, pG, pH, pL checked so far).
- [ ] Re-run the DCC/EOP association scripts against the combined 317-unit
  typing table (`plasmid_final_types_with_hypothetical_v2.tsv`) rather than
  the original 210-unit table, now that hypothetical plasmids are typed.
