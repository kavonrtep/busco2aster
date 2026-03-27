rule infer_gene_tree:
    input:
        alignment=f"{ALIGNMENT_DIR}/{{locus_id}}.aln.faa",
    output:
        command="results/gene_trees/per_locus/{locus_id}/command.sh",
        report="results/gene_trees/per_locus/{locus_id}/{locus_id}.iqtree",
        log="results/gene_trees/per_locus/{locus_id}/{locus_id}.log",
        treefile="results/gene_trees/per_locus/{locus_id}/{locus_id}.treefile",
        completion="results/gene_trees/per_locus/{locus_id}/run.complete",
    params:
        executable=config["iqtree_executable"],
        prefix=lambda wildcards: gene_tree_output_paths(wildcards.locus_id)["prefix"],
        model=config["iqtree_model"],
        support_mode=config["iqtree_support_mode"],
        seed=int(config["iqtree_seed"]),
        ufboot_replicates=int(config["iqtree_ufboot_replicates"]),
    threads:
        get_thread_count("iqtree")
    run:
        import subprocess
        from pathlib import Path

        from scripts.gene_trees import build_iqtree_command, write_command_script

        command = build_iqtree_command(
            executable=params.executable,
            alignment_path=input.alignment,
            prefix=params.prefix,
            threads=threads,
            model=params.model,
            support_mode=params.support_mode,
            seed=params.seed,
            ufboot_replicates=params.ufboot_replicates,
        )
        write_command_script(Path(output.command), command)
        subprocess.run(command, check=True)
        Path(output.completion).touch()


rule aggregate_gene_trees:
    input:
        retained=RETAINED_LOCI_TABLE,
        treefiles=retained_gene_tree_treefiles,
        reports=retained_gene_tree_reports,
    output:
        manifest=GENE_TREE_MANIFEST,
        aggregate=GENE_TREE_AGGREGATE,
    params:
        support_mode=config["iqtree_support_mode"],
    shell:
        (
            "python3 -m scripts.aggregate_gene_trees "
            "--retained {input.retained:q} "
            "--manifest {output.manifest:q} "
            "--output {output.aggregate:q} "
            "--support-mode {params.support_mode:q}"
        )
