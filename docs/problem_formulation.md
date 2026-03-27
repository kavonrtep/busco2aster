Goal
Build a reproducible phylogenomics workflow that takes genome assemblies as input and returns an unrooted coalescent-aware species tree, using BUSCO-defined ortholog markers, per-locus ML gene trees from IQ-TREE 2, and species-tree inference with ASTRAL-III. The workflow must also emit enough QC to explain which taxa and loci were kept, filtered, or flagged.

1) What the pipeline is actually solving

You are not just “building a tree from BUSCO genes.” You are solving a more specific problem:

identify a comparable set of near-universal ortholog markers across assemblies,
extract one homologous sequence per retained taxon and locus,
infer one gene tree per locus,
account for gene-tree discordance under the multispecies coalescent when estimating the species tree.

That framing fits the tools you picked. BUSCO lineage datasets are built from marker genes that are present in at least 90% of species in a lineage and single-copy in 90% of those species; BUSCO genome mode works directly from nucleotide assemblies and classifies hits as complete single-copy, duplicated, fragmented, or missing. ASTRAL-III is explicitly a summary method that estimates an unrooted species tree from a set of unrooted gene trees under the MSC/quartet framework.

2) Design decisions to settle now
A. What counts as a valid marker set

Use one BUSCO lineage dataset shared by all taxa in the study. Do not let different samples choose different lineages if the goal is cross-sample ortholog comparability. BUSCO is lineage-specific, and BUSCO v6 is the current stable release.

My recommendation: pick the most specific lineage that still covers every species in your analysis. That is an implementation choice, not a BUSCO rule, but it usually gives a better marker set than using a very broad lineage.

B. Genome-mode sequence source

This is the first place your pipeline can quietly go wrong. In BUSCO genome mode, eukaryotic runs now default to Miniprot, which BUSCO describes as a gene mapper rather than a gene predictor. BUSCO also warns that, for divergent assemblies, Miniprot may underperform and some reported gene predictions may contain internal stop codons and may not be suitable for downstream analysis.

So for implementation, I would treat BUSCO-derived sequences as candidates, not unquestioned truth. Add a sanity-check layer for stop codons, excessive truncation, suspicious length outliers, and abnormal gap content.

C. Missing data policy

Do not hard-code “locus must be present in all taxa.” That sounds clean, but it is often the wrong choice. ASTRAL can use gene trees with missing taxa, and work summarized in the ASTRAL documentation cites evidence that filtering genes purely because of missing data is often neutral or harmful to species-tree accuracy.

My default would be an occupancy threshold such as 70–90% taxa per locus, not 100%.

D. Duplication policy

For the first version, be strict: retain only loci where included taxa are complete and single-copy. BUSCO explicitly separates single_copy_busco_sequences from multi_copy_busco_sequences, so the data structure already supports this.

But plan for failure here. If many taxa are duplicated because of WGD, assembly redundancy, haplotig retention, or real paralogy, classical ASTRAL-III is no longer the natural endpoint. The ASTRAL project now points users to ASTRAL-Pro for multi-copy genes, and the newer ASTER family also includes Weighted ASTRAL.

E. Gene-tree uncertainty

Do not send raw, fully resolved, weak gene trees straight into ASTRAL and assume the species tree will fix it. ASTRAL’s own tutorial recommends contracting very low-support branches in gene trees because that can improve accuracy; their example threshold is branches below 10% bootstrap support. IQ-TREE also warns that UFBoot support is not interpreted the same way as standard bootstrap, and recommends at least 1000 UFBoot replicates; -bnni is advised when model violations may inflate support.

3) Workflow I would implement
Stage 0 — Inputs and metadata

Required inputs:

avaialable test_data - in @test_data directory 

genome assembly FASTA per sample
sample sheet: sample_id, species_id, optional group, optional outgroup
chosen BUSCO lineage
global settings: occupancy threshold, sequence type, support-collapse threshold

Strongly recommended:

one genome per species for v1
sanitized taxon names from the start

This matters because IQ-TREE normalizes problematic characters in sequence names, and ASTRAL is picky about taxon names as well.

Stage 1 — BUSCO per assembly

Run BUSCO in genome mode on every assembly using the same lineage. Collect:

short_summary.*
full_table.tsv
busco_sequences/single_copy_busco_sequences
duplicated/fragmented/missing counts

BUSCO’s output structure already exposes these pieces, including single_copy_busco_sequences, multi_copy_busco_sequences, and full_table.tsv.

Stage 2 — Build the locus-by-taxon matrix

For every BUSCO ID, build a matrix with one row per locus and one column per taxon, storing:

status: single / duplicated / fragmented / missing
sequence length
translated/protein length if relevant
stop-codon / invalid-character flags
source file path

This matrix is the core object in your pipeline. Everything downstream should derive from it.

Stage 3 — Locus filtering

Suggested v1 filtering logic:

keep only BUSCO loci with taxon occupancy ≥ threshold
exclude any taxon-locus cell not classified as complete single-copy
optionally drop loci with extreme length variance
optionally drop individual sequences with clearly fragmentary or pathological signal

This is where I would be selective, but not aggressive. Evidence cited by the ASTRAL authors argues that removing fragmentary sequences is useful, whereas removing loci just because of missing data is often not.

Stage 4 — Per-locus sequence preparation

Create one FASTA per BUSCO locus from retained taxa only.

For v1, I would use protein sequences as the default analysis layer. That is my implementation recommendation, not a rule from the docs. It is usually the safer starting point when the upstream sequences come from genome-based ortholog extraction and divergence may be nontrivial.

Add optional modes later:

CDS nucleotide alignment
codon-aware back-translation
separate shallow/deep divergence strategies
Stage 5 — Per-locus alignment and optional trimming

Align each locus independently. Make trimming optional and conservative. There is no universal win from aggressive trimming; newer work like ClipKIT was proposed specifically because conventional trimming can worsen phylogenetic inference, whereas retaining phylogenetically informative sites can be more robust.

So the right design is:

align every locus independently
keep raw alignment
optionally generate a trimmed alignment
record which one was actually used for tree inference
Stage 6 — Gene trees with IQ-TREE 2

Infer one ML tree per locus with per-locus model selection. IQ-TREE 2 added direct support for directories of single-locus alignments, and IQ-TREE’s concordance-factor workflow documents exactly this pattern: infer separate locus trees from a folder and collect them in one loci.treefile.

For the pipeline spec, each locus should produce at least:

best-fit model
ML tree
support metrics
log/report file

I would make support estimation part of the default run, because weak-branch contraction is easier if support is already there.

Stage 7 — Collapse weak gene-tree branches

Before species-tree inference, convert obviously unreliable gene-tree edges into polytomies. This is not optional in spirit, even if you keep the threshold configurable. ASTRAL can accept unresolved gene trees, and its authors explicitly recommend removing very low-support branches.

Stage 8 — Species tree with ASTRAL-III

Feed the collapsed gene trees to ASTRAL-III. ASTRAL estimates an unrooted species tree from unrooted gene trees; the input trees may contain missing taxa and polytomies. The output should also be treated as unrooted.

For implementation, produce two ASTRAL outputs:

the inferred species tree
a scored/annotated version using -q / richer branch annotations

ASTRAL reports branch length in coalescent units and local posterior support, and with richer annotations it can also report quartet-frequency quantities and effective number of genes (EN) per branch.

Stage 9 — Post hoc diagnostics

This is where the pipeline becomes useful instead of just runnable.

I would calculate:

ASTRAL normalized quartet score
local posterior support
EN per branch
IQ-TREE gCF and sCF on the final species tree

IQ-TREE explicitly allows gCF/sCF calculation on any reference tree, including one inferred by a coalescent method, and it documents the standard workflow for this.

That gives you three layers of evidence:

branch support from ASTRAL
how many genes are informative for a branch (EN)
how much the loci and sites actually agree with that branch (gCF/sCF)
4) What the pipeline should output

I would define these as required deliverables:

busco_summary.tsv — one row per sample
locus_taxon_matrix.tsv — one row per BUSCO locus, one column block per taxon
retained_loci.tsv — why each locus passed or failed
per_locus_alignments/
per_locus_trees/
gene_trees.raw.tre
gene_trees.collapsed.tre
species_tree.astral.tre
species_tree.astral.annotated.tre
species_tree.astral.log
species_tree.gcf_scf.tre
pipeline_report.html or report.md

That report should explicitly summarize:

number of input taxa
BUSCO completeness per taxon
number of candidate loci
number retained after filtering
occupancy distribution
duplicated-locus rate
branch-support summary
discordance summary
5) What I would not do in v1

I would avoid these in the first implementation:

mixed BUSCO lineages across taxa
automatic switching between protein and nucleotide modes
hidden heuristics for duplicated BUSCO rescue
per-clade custom thresholds
rooting inside the ASTRAL stage
trying to solve paralogy with ad hoc rules

Keep v1 strict and inspectable.

6) The v1 spec I’d actually ship

If I had to freeze a first-pass design now, it would be:

Input: one genome assembly per species
BUSCO: one shared lineage, genome mode
Loci: complete single-copy only
Occupancy threshold: 80%
Sequence type: amino acid
Alignment: per locus, optional conservative trimming
Gene trees: IQ-TREE 2 with per-locus model selection and support
Uncertainty handling: collapse weak branches before ASTRAL
Species tree: ASTRAL-III
Diagnostics: ASTRAL support/quartet score + IQ-TREE gCF/sCF
Output: unrooted species tree plus full locus/QC audit trail
7) One forward-looking change worth planning for

Even if you start with ASTRAL-III, I would design the species-tree stage as a backend interface, not hardwired logic. The ASTRAL/ASTER ecosystem now includes ASTRAL-Pro for multi-copy gene families and Weighted ASTRAL for gene trees with branch lengths/supports. That makes it easy to keep your strict BUSCO + ASTRAL-III workflow now, but not trap yourself later if duplication or scale becomes the real bottleneck.

If you want, the next step should be turning this into a concrete pipeline spec with:

directory structure,
exact input/output filenames,
per-rule filtering logic,
and a Snakemake  implementatio
conda environments - if possible - this need to be checked - IQ tree 2(or 3) is probalbly available as binaries, 
but ASTRAL-III / or better ASTER ?  - https://github.com/chaoszhang/ASTER



