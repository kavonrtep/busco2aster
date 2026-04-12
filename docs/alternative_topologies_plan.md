# Alternative Topology Hypothesis Testing — Problem Formulation and Implementation Plan

## 1. Problem Statement

The current `busco2aster` pipeline produces a single species tree — the one
that maximizes the wASTRAL quartet score given the input gene trees. The
concordance analysis (gCF, sCF, wASTRAL quartet frequencies) already quantifies
*per-branch* support, but the pipeline does not answer a natural follow-up
question: **which alternative species tree topologies are also consistent with
the data, and how much worse are they compared to the best tree?**

This matters in practice because:

- Many phylogenomic datasets contain branches where the best resolution is only
  marginally better than the alternative. Reporting only the best tree hides
  this uncertainty.
- Reviewers and collaborators increasingly expect formal topology testing (AU
  test, SH test) alongside the point-estimate species tree.
- Identifying contested branches and their alternatives clarifies whether
  discordance is driven by ILS, hybridization, or insufficient data — each
  calling for different biological interpretation.

The goal of this extension is to add a principled, automated workflow that:

1. Identifies branches in the species tree where the topology is uncertain.
2. Generates a small set of plausible alternative species tree topologies.
3. Formally tests whether the alternatives can be statistically rejected.
4. Reports the results in the existing HTML/Markdown report framework.


## 2. Background: Sources of Topological Uncertainty

There are three distinct layers where uncertainty enters the pipeline. The
extension must address layers 2 and 3; layer 1 is already partially handled
by the existing branch support (aBayes) on gene trees.

### 2.1 Gene tree estimation error (layer 1)

Each gene tree is a point estimate from IQ-TREE. Finite alignment length and
model misspecification mean the inferred gene tree may differ from the true
gene tree. The `abayes` branch support values already quantify this per-branch
in each gene tree, and wASTRAL uses these support values as weights — gene tree
quartets from poorly supported branches contribute less to the species tree
score. No additional work is needed here.

### 2.2 Gene tree discordance due to ILS (layer 2)

Under the multispecies coalescent (MSC), different genomic regions have
genuinely different evolutionary histories due to incomplete lineage sorting
(ILS). This is not error — it is expected biology. wASTRAL accounts for this
by finding the species tree topology that best explains the distribution of
quartet topologies across gene trees.

The key insight is that for each internal branch in the species tree, the four
subtrees hanging off that branch define three possible unrooted quartet
resolutions. Under the MSC, the species tree topology always produces the most
frequent quartet — but the margin can be slim when internal branches are short
(rapid radiations, recent divergences). The `-u 2` annotation mode in wASTRAL
reports the quartet support and local posterior probability for all three
resolutions of every branch, directly quantifying this uncertainty.

### 2.3 Species tree estimation uncertainty (layer 3)

Even given the gene trees, the species tree inference is a statistical
estimation. Two topologies might have nearly identical wASTRAL scores. The
AU test on a concatenated supermatrix provides a complementary,
likelihood-based assessment of whether two candidate topologies can be
distinguished given the sequence data.

The AU test (Shimodaira, 2002) compares log-likelihoods of candidate trees
evaluated on the same alignment using a multiscale bootstrap resampling
procedure. Trees with AU p-value >= 0.05 belong to the confidence set and
cannot be statistically rejected.


## 3. Approach Overview

The extension adds three sequential analysis stages after the existing species
tree inference:

```
existing pipeline
─────────────────
gene trees ──► wASTRAL species tree ──► concordance (gCF, sCF)

new extension
─────────────
         ├──► [A] Contested branch identification (wASTRAL -u 2)
         ├──► [B] Alternative topology generation (NNI at contested branches)
         └──► [C] AU topology test (IQ-TREE on concatenated supermatrix)
                    │
                    ▼
              extended report
```

### Stage A: Contested branch identification

Run wASTRAL with `-u 2` to annotate every internal branch with:

- `q1`, `q2`, `q3` — quartet support for the three resolutions (fraction of
  weighted gene tree quartets supporting each topology)
- `pp1`, `pp2`, `pp3` — local posterior probabilities for the three resolutions
  (Bayesian conversion of quartet scores; sum to 1.0)

A branch is flagged as **contested** if `pp1 < contested_branch_threshold`
(configurable, default `0.95`).

Output: a TSV table listing every internal branch with its quartet support
values, posterior probabilities, and a boolean `contested` flag.

### Stage B: Alternative topology generation

For each contested branch, the alternative topology is already specified by the
`q2`/`q3` annotations — it is the NNI swap that reconnects the four subtrees
(L, R, S, O) as LS|RO or LO|RS instead of the current LR|SO.

Alternative species trees are generated by applying NNI swaps at contested
branches, one at a time and in combinations:

- With k contested branches, generate up to 2^k - 1 alternative trees (each
  branch can be either in its original or swapped state, excluding the original
  tree itself).
- In practice, k is typically 1–3 for well-sampled phylogenomic datasets,
  yielding 1–7 alternatives.
- If k > `max_contested_branches` (configurable, default `4`), only the top-k
  most contested branches (sorted by ascending pp1) are used, and a warning is
  emitted.

Additionally, users can supply arbitrary hypothesis trees via a config option
`hypothesis_trees` pointing to a Newick file. These are appended to the
candidate set without modification.

Output: a multi-tree Newick file containing the wASTRAL tree (first) followed
by all generated alternatives and user-supplied hypotheses.

### Stage C: AU topology test

The candidate trees from Stage B are tested against the sequence data using the
IQ-TREE approximately unbiased (AU) test. This requires:

1. **Supermatrix construction**: concatenate the per-locus alignments into a
   single alignment with a partition file defining locus boundaries. Missing
   taxa are filled with gap characters. The partition file assigns each locus
   its own substitution model (or uses the same model selected during gene tree
   inference).

2. **AU test execution**: IQ-TREE evaluates the log-likelihood of each
   candidate tree on the supermatrix and performs the AU test via RELL
   resampling.

   ```
   iqtree3 -s supermatrix.phy \
           -p partitions.nex \
           --trees candidates.tre \
           --test-au \
           -n 0 \
           --test 10000 \
           --prefix au_test
   ```

   The `-n 0` flag tells IQ-TREE not to search for a new ML tree but only to
   evaluate the provided trees (optimizing branch lengths per tree).

3. **Result parsing**: extract the per-tree log-likelihood, delta-L, and AU
   p-value from the `.iqtree` output file. Trees with p-AU >= 0.05 are in the
   confidence set.

Output: a TSV with one row per candidate tree, columns for log-likelihood,
delta-L, and AU p-value, plus a boolean `in_confidence_set` flag.


## 4. Configuration

New config keys in `config/config.yaml`:

```yaml
# ── Alternative topology testing ─────────────────────────────────────────
# Enable/disable the entire extension (default: true when pipeline v2+)
run_topology_tests: true

# pp1 threshold below which a branch is flagged as contested
contested_branch_threshold: 0.95

# Maximum number of contested branches to use for alternative generation
# (branches are ranked by ascending pp1; excess branches are ignored)
max_contested_branches: 4

# Optional: user-supplied hypothesis trees (Newick file, one tree per line)
# These are added to the candidate set alongside generated alternatives.
hypothesis_trees: null     # e.g. "config/hypothesis_trees.tre"

# Number of RELL replicates for the AU test (default: 10000)
au_test_replicates: 10000

# Substitution model for AU test supermatrix partitions.
# "from_gene_trees" re-uses the model selected during gene tree inference.
# Otherwise, specify a model string (e.g. "MFP" for automatic selection).
au_test_model: "from_gene_trees"
```


## 5. Implementation Plan

### 5.1 New output paths

```
results/
├── topology_tests/
│   ├── branch_quartet_support.tsv       # Stage A output
│   ├── contested_branches.tsv           # Stage A: filtered subset
│   ├── candidate_trees.tre              # Stage B: multi-tree Newick
│   ├── candidate_trees_manifest.tsv     # Stage B: tree index → description
│   ├── supermatrix.phy                  # Stage C input: concatenated alignment
│   ├── supermatrix_partitions.nex       # Stage C input: partition definitions
│   ├── au_test.iqtree                   # Stage C: IQ-TREE report
│   ├── au_test_results.tsv             # Stage C: parsed AU test results
│   └── topology_tests.complete          # sentinel
└── report/
    └── data/
        ├── branch_quartet_support.tsv   # copy for report
        ├── candidate_trees_manifest.tsv # copy for report
        └── au_test_results.tsv          # copy for report
```

### 5.2 New scripts

**`scripts/quartet_support.py`** — Stage A

- Input: wASTRAL annotated tree (Newick with `-u 2` bracket annotations).
- Parse the tree using `ete3` or `dendropy`, extract per-branch `q1/q2/q3`
  and `pp1/pp2/pp3` from the annotation strings.
- Write `branch_quartet_support.tsv` with columns: `branch_id`, `taxa_left`,
  `taxa_right`, `taxa_sister`, `taxa_outgroup`, `q1`, `q2`, `q3`, `pp1`,
  `pp2`, `pp3`, `contested`.
- Write `contested_branches.tsv` (filtered rows where `pp1 < threshold`).

**`scripts/generate_alternatives.py`** — Stage B

- Input: wASTRAL tree, contested branches TSV, optional user hypothesis trees.
- For each contested branch, perform NNI swap: detach the subtree
  corresponding to the best alternative (q2 or q3, whichever is higher) and
  reattach in the alternative position.
- Generate all combinations of swaps up to `max_contested_branches`.
- Write `candidate_trees.tre` (multi-tree Newick, original tree first).
- Write `candidate_trees_manifest.tsv` with columns: `tree_index`,
  `description` (e.g., "wASTRAL best", "branch_3 swapped to alt1",
  "branch_3+branch_7 swapped"), `source` (generated / user_hypothesis).

**`scripts/build_supermatrix.py`** — Stage C input preparation

- Input: per-locus alignment directory, retained loci TSV, gene tree model
  info (from IQ-TREE `.iqtree` or `.best_model.nex` files).
- Read each alignment, pad taxa missing from that locus with all-gap
  sequences, concatenate.
- Write `supermatrix.phy` in relaxed PHYLIP format.
- Write `supermatrix_partitions.nex` as a NEXUS charsets block with per-locus
  model assignments.
- Taxon ordering is consistent with the species tree leaf labels.

**`scripts/parse_au_test.py`** — Stage C output parsing

- Input: IQ-TREE `.iqtree` report file from the AU test run.
- Parse the topology test table (log-likelihood, delta-L, p-AU, KH, SH
  columns).
- Join with `candidate_trees_manifest.tsv` to attach tree descriptions.
- Write `au_test_results.tsv`.

### 5.3 New Snakemake rules

**`workflow/rules/topology_tests.smk`** — single new rule file containing:

| Rule | Depends on | Produces | Notes |
|------|-----------|----------|-------|
| `wastral_annotated` | gene trees (existing) | annotated wASTRAL tree | Re-runs wASTRAL with `-u 2`. Can reuse the existing species tree topology and just add annotation, or re-run from scratch (cheap). |
| `identify_contested_branches` | annotated tree | `branch_quartet_support.tsv`, `contested_branches.tsv` | Calls `scripts/quartet_support.py`. Pure Python, no external tools. |
| `generate_alternative_trees` | wASTRAL tree, contested branches, config | `candidate_trees.tre`, `candidate_trees_manifest.tsv` | Calls `scripts/generate_alternatives.py`. Pure Python. |
| `build_supermatrix` | per-locus alignments, retained loci, gene tree models | `supermatrix.phy`, `supermatrix_partitions.nex` | Calls `scripts/build_supermatrix.py`. I/O-bound, no external tools. |
| `au_topology_test` | supermatrix, partitions, candidate trees | `au_test.iqtree`, `au_test_results.tsv` | Runs IQ-TREE with `--trees --test-au -n 0`. Uses `IQTREE_EXECUTABLE`. Most expensive step but still fast (branch-length optimization only, no tree search). |
| `topology_tests_complete` | all above | `topology_tests.complete` | Sentinel rule. |

### 5.4 Modifications to existing code

| File | Change |
|------|--------|
| `Snakefile` | Add `include: "workflow/rules/topology_tests.smk"`. Add topology test outputs to `rule all` inputs (gated by `config["run_topology_tests"]`). Add new output path constants. |
| `config/config.yaml` | Add new configuration keys (section 4). |
| `config/config.template.yaml` | Add new keys with default values and comments. |
| `workflow/rules/report.smk` | Add topology test TSVs as inputs to the report data collection rule. |
| `workflow/rules/visual_report.smk` | Add new report section rendering the topology test results (contested branches table, candidate tree descriptions, AU test results table, confidence set summary). |
| `reports/` (Quarto/R templates) | Add a new section "Topological uncertainty" with: (1) table of per-branch quartet posteriors with contested branches highlighted, (2) depiction of alternative trees showing which branches differ, (3) AU test results table with confidence set annotation. |

### 5.5 Dependencies

No new external tool dependencies — the extension uses wASTRAL (already
installed) and IQ-TREE (already installed). The Python scripts use only
standard library plus `ete3` or `dendropy` for tree manipulation, which should
be added to the pipeline's Python environment if not already present.

The Singularity container definition (`busco2aster.def`) needs to include the
tree manipulation library. If `ete3` is too heavy (it pulls in PyQt), `dendropy`
is the lighter alternative and sufficient for NNI operations.

### 5.6 Computational cost

The extension adds negligible cost to the overall pipeline:

- **wASTRAL with `-u 2`**: same runtime as the existing wASTRAL call (seconds
  to minutes, depending on taxon count).
- **Supermatrix construction**: I/O-bound, typically seconds.
- **AU test**: IQ-TREE evaluates each candidate tree on the supermatrix by
  optimizing branch lengths only (no tree search). For a typical BUSCO-scale
  dataset (~500 loci, ~20 taxa, ~5 candidate trees), this takes minutes, not
  hours. The RELL resampling (10,000 replicates) adds modest overhead.

The total added wall-clock time is on the order of the existing concordance
analysis — a small fraction of the gene tree inference step that dominates
the pipeline runtime.


## 6. Interpretation Guide (for report documentation)

The report section should include a brief interpretation guide for users:

- **Quartet posteriors (pp1, pp2, pp3)**: These are coalescent-model-aware.
  They account for ILS and answer: "Given the gene tree discordance pattern,
  how probable is each resolution of this branch under the multispecies
  coalescent?" A branch with pp1 > 0.95 is well-resolved. A branch with
  pp1 ~ 0.5–0.7 indicates genuine topological ambiguity.

- **AU test p-values**: These are concatenation-based. They answer: "Can we
  reject this whole-tree topology based on the total sequence evidence?" Trees
  with p-AU >= 0.05 cannot be rejected and belong to the confidence set.

- **Agreement between the two**: When a contested branch (low pp1) corresponds
  to an alternative tree that is also in the AU confidence set (p-AU >= 0.05),
  the topological uncertainty is robust — both the coalescent and concatenation
  perspectives agree that the data cannot resolve this branch. This is the
  strongest signal of genuine ambiguity.

- **Disagreement**: If pp1 is low but the AU test rejects the alternative, the
  gene tree discordance may be due to systematic error (e.g., long-branch
  attraction, model misspecification) rather than ILS. Conversely, if pp1 is
  high but the AU test does not reject an alternative, ILS is absorbing
  discordance that the concatenation approach cannot model.


## 7. Testing Strategy

- **Unit tests for new scripts**: quartet support parser, NNI swap correctness
  (verify that swapping a branch and swapping back recovers the original tree),
  supermatrix construction (verify alignment lengths sum correctly, gap padding
  is correct).
- **Integration test with test dataset**: run the full extension on the
  existing test data. Verify that all output files are produced, that the
  wASTRAL tree appears as the first tree in `candidate_trees.tre`, and that the
  AU test `.iqtree` file contains the expected number of trees.
- **Edge cases**: zero contested branches (no alternatives generated, AU test
  runs on the wASTRAL tree alone as a degenerate case), all branches contested
  (cap at `max_contested_branches`), single-taxon loci in the supermatrix
  (should be excluded or gap-padded).


## 8. References

- Shimodaira, H. (2002). An approximately unbiased test of phylogenetic tree
  selection. Systematic Biology, 51(3), 492–508.
- Zhang, C., Rabiee, M., Sayyari, E., & Mirarab, S. (2018). ASTRAL-III:
  polynomial time species tree reconstruction from partially resolved gene
  trees. BMC Bioinformatics, 19(S6), 153.
- Zhang, C., & Mirarab, S. (2022). Weighting by gene tree uncertainty improves
  accuracy of quartet-based species trees. Molecular Biology and Evolution,
  39(12), msac215.
- Zhang, C., Nielsen, R., & Mirarab, S. (2025). ASTER: A Package for
  Large-scale Phylogenomic Reconstructions. Molecular Biology and Evolution,
  42(8), msaf172.
- Minh, B. Q., Schmidt, H. A., Chernomor, O., et al. (2020). IQ-TREE 2: New
  models and efficient methods for phylogenetic inference in the genomic era.
  Molecular Biology and Evolution, 37(5), 1530–1534.
- Stenz, N. W. M., et al. (2023). Likelihood-based tests of species tree
  hypotheses. Systematic Biology, 72(5).
