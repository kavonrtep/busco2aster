rule wastral_u2_annotate:
    """Re-run wASTRAL with -u 2 to annotate every branch with quartet posteriors (pp1/pp2/pp3)."""
    input:
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        gene_trees=WASTRAL_GENE_TREE_INPUT,
    output:
        command=TOPOLOGY_TESTS_OUTPUTS["wastral_u2_command"],
        tree=TOPOLOGY_TESTS_OUTPUTS["wastral_u2_tree"],
        log=TOPOLOGY_TESTS_OUTPUTS["wastral_u2_log"],
    params:
        executable=WASTRAL_EXECUTABLE,
        support_mode=IQTREE_SUPPORT_MODE,
        annotation_mode=2,
    threads:
        get_thread_count("species_tree")
    run:
        import shutil
        import subprocess
        import tempfile
        from pathlib import Path

        from scripts.species_tree import build_aster_command, write_species_tree_command_script

        work_root = Path("work") / "topology_tests"
        work_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="wastral_u2.", dir=work_root.as_posix()) as tmpdir:
            tmp_path = Path(tmpdir).resolve()
            tree_output = tmp_path / "annotated.u2.tre"
            command = build_aster_command(
                backend="wastral",
                executable=Path(params.executable).resolve().as_posix(),
                input_path=Path(input.gene_trees).resolve().as_posix(),
                output_path=tree_output.as_posix(),
                threads=threads,
                support_mode=params.support_mode,
                score_constraint_tree=True,
                constraint_path=Path(input.species_tree).resolve().as_posix(),
                annotation_mode=params.annotation_mode,
            )
            write_species_tree_command_script(
                Path(output.command),
                command,
                Path(output.log),
                cwd=tmp_path,
            )
            with Path(output.log).open("w", encoding="utf-8") as log_handle:
                subprocess.run(command, check=True, cwd=tmp_path, stderr=log_handle)
            if not tree_output.is_file():
                raise FileNotFoundError(
                    f"wASTRAL -u 2 annotation did not write expected tree: {tree_output}"
                )
            Path(output.tree).parent.mkdir(parents=True, exist_ok=True)
            shutil.move(tree_output.as_posix(), output.tree)


rule identify_contested_branches:
    """Parse wASTRAL -u 2 tree and write per-branch quartet posteriors + contested subset."""
    input:
        u2_tree=TOPOLOGY_TESTS_OUTPUTS["wastral_u2_tree"],
    output:
        branch_quartet_support=TOPOLOGY_TESTS_OUTPUTS["branch_quartet_support"],
        contested_branches=TOPOLOGY_TESTS_OUTPUTS["contested_branches"],
    params:
        threshold=CONTESTED_BRANCH_THRESHOLD,
    run:
        from pathlib import Path

        from scripts.quartet_support import (
            filter_contested_branches,
            parse_wastral_u2_tree,
            write_branch_quartet_support,
            write_contested_branches,
        )

        rows = parse_wastral_u2_tree(
            Path(input.u2_tree),
            contested_threshold=params.threshold,
        )
        contested = filter_contested_branches(rows)
        write_branch_quartet_support(rows, Path(output.branch_quartet_support))
        write_contested_branches(contested, Path(output.contested_branches))


rule generate_alternative_trees:
    """Generate NNI-based alternative species tree topologies at contested branches."""
    input:
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        contested_branches=TOPOLOGY_TESTS_OUTPUTS["contested_branches"],
    output:
        candidate_trees=TOPOLOGY_TESTS_OUTPUTS["candidate_trees"],
        candidate_manifest=TOPOLOGY_TESTS_OUTPUTS["candidate_manifest"],
    params:
        max_contested_branches=MAX_CONTESTED_BRANCHES,
        hypothesis_trees=HYPOTHESIS_TREES,
    run:
        import csv
        from pathlib import Path

        from scripts.generate_alternatives import (
            generate_candidate_trees,
            write_candidate_manifest,
            write_candidate_trees,
        )

        with open(input.contested_branches, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            contested_rows = [dict(row) for row in reader]

        hypothesis_path = (
            Path(params.hypothesis_trees)
            if params.hypothesis_trees is not None
            else None
        )

        trees, descriptions, sources = generate_candidate_trees(
            species_tree_path=Path(input.species_tree),
            contested_rows=contested_rows,
            max_contested_branches=params.max_contested_branches,
            hypothesis_path=hypothesis_path,
        )
        write_candidate_trees(trees, Path(output.candidate_trees))
        write_candidate_manifest(descriptions, sources, Path(output.candidate_manifest))


rule build_supermatrix:
    """Concatenate per-locus alignments into a PHYLIP supermatrix for AU testing."""
    input:
        alignments_complete=ALIGNMENTS_COMPLETE,
        retained_loci=RETAINED_LOCI_TABLE,
        gene_tree_manifest=GENE_TREE_MANIFEST,
    output:
        supermatrix=TOPOLOGY_TESTS_OUTPUTS["supermatrix"],
        partitions=TOPOLOGY_TESTS_OUTPUTS["partitions"],
    params:
        alignment_dir=ALIGNMENT_DIR,
        alignment_suffix=f".{ALIGNMENT_SUFFIX}",
        au_test_model=AU_TEST_MODEL,
    run:
        from pathlib import Path

        from scripts.build_supermatrix import build_supermatrix

        build_supermatrix(
            alignment_dir=Path(params.alignment_dir),
            retained_loci_path=Path(input.retained_loci),
            gene_tree_manifest_path=Path(input.gene_tree_manifest),
            output_phy=Path(output.supermatrix),
            output_nex=Path(output.partitions),
            au_test_model=params.au_test_model,
            alignment_suffix=params.alignment_suffix,
        )


rule au_topology_test:
    """Run IQ-TREE AU test on candidate trees using the concatenated supermatrix."""
    input:
        supermatrix=TOPOLOGY_TESTS_OUTPUTS["supermatrix"],
        partitions=TOPOLOGY_TESTS_OUTPUTS["partitions"],
        candidate_trees=TOPOLOGY_TESTS_OUTPUTS["candidate_trees"],
        candidate_manifest=TOPOLOGY_TESTS_OUTPUTS["candidate_manifest"],
    output:
        command=TOPOLOGY_TESTS_OUTPUTS["au_command"],
        au_iqtree=TOPOLOGY_TESTS_OUTPUTS["au_iqtree"],
        au_results=TOPOLOGY_TESTS_OUTPUTS["au_results"],
        au_log=TOPOLOGY_TESTS_OUTPUTS["au_log"],
    params:
        executable=IQTREE_EXECUTABLE,
        prefix=TOPOLOGY_TESTS_OUTPUTS["au_prefix"],
        replicates=AU_TEST_REPLICATES,
    threads:
        get_thread_count("concordance")
    run:
        import subprocess
        from pathlib import Path

        from scripts.parse_au_test import join_with_manifest, parse_au_test_iqtree, write_au_results
        from scripts.topology_tests import build_au_test_command, write_topology_tests_command_script

        command = build_au_test_command(
            executable=params.executable,
            supermatrix=input.supermatrix,
            partitions=input.partitions,
            candidate_trees=input.candidate_trees,
            prefix=params.prefix,
            threads=threads,
            replicates=params.replicates,
        )
        write_topology_tests_command_script(
            Path(output.command),
            command,
            Path(output.au_log),
        )
        with Path(output.au_log).open("w", encoding="utf-8") as log_handle:
            subprocess.run(command, check=True, stderr=log_handle)

        au_rows = parse_au_test_iqtree(Path(output.au_iqtree))
        joined = join_with_manifest(au_rows, Path(input.candidate_manifest))
        write_au_results(joined, Path(output.au_results))


rule topology_tests_complete:
    """Sentinel rule: all topology test outputs are ready."""
    input:
        branch_quartet_support=TOPOLOGY_TESTS_OUTPUTS["branch_quartet_support"],
        contested_branches=TOPOLOGY_TESTS_OUTPUTS["contested_branches"],
        candidate_trees=TOPOLOGY_TESTS_OUTPUTS["candidate_trees"],
        candidate_manifest=TOPOLOGY_TESTS_OUTPUTS["candidate_manifest"],
        supermatrix=TOPOLOGY_TESTS_OUTPUTS["supermatrix"],
        partitions=TOPOLOGY_TESTS_OUTPUTS["partitions"],
        au_iqtree=TOPOLOGY_TESTS_OUTPUTS["au_iqtree"],
        au_results=TOPOLOGY_TESTS_OUTPUTS["au_results"],
    output:
        touch(TOPOLOGY_TESTS_OUTPUTS["completion"])
