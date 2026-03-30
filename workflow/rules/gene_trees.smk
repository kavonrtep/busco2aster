rule infer_gene_trees:
    input:
        alignments_complete=ALIGNMENTS_COMPLETE,
    output:
        command=GENE_TREE_DIRECTORY_OUTPUTS["command"],
        report=GENE_TREE_DIRECTORY_OUTPUTS["report"],
        log=GENE_TREE_DIRECTORY_OUTPUTS["log"],
        treefile=GENE_TREE_DIRECTORY_OUTPUTS["treefile"],
        best_model_nex=GENE_TREE_DIRECTORY_OUTPUTS["best_model_nex"],
        best_scheme=GENE_TREE_DIRECTORY_OUTPUTS["best_scheme"],
        best_scheme_nex=GENE_TREE_DIRECTORY_OUTPUTS["best_scheme_nex"],
        model_gz=GENE_TREE_DIRECTORY_OUTPUTS["model_gz"],
        checkpoint=GENE_TREE_DIRECTORY_OUTPUTS["checkpoint"],
        parstree=GENE_TREE_DIRECTORY_OUTPUTS["parstree"],
    params:
        executable=IQTREE_EXECUTABLE,
        prefix=GENE_TREE_DIRECTORY_OUTPUTS["prefix"],
        alignment_dir=ALIGNMENT_DIR,
        model=IQTREE_MODEL,
        support_mode=IQTREE_SUPPORT_MODE,
        seed=IQTREE_SEED,
        ufboot_replicates=IQTREE_UFBOOT_REPLICATES,
    threads:
        get_thread_count("iqtree")
    run:
        import subprocess
        from pathlib import Path

        from scripts.gene_trees import build_iqtree_directory_command, write_command_script

        command = build_iqtree_directory_command(
            executable=params.executable,
            alignment_dir=params.alignment_dir,
            prefix=params.prefix,
            threads=threads,
            model=params.model,
            support_mode=params.support_mode,
            seed=params.seed,
            ufboot_replicates=params.ufboot_replicates,
        )
        write_command_script(Path(output.command), command)
        with Path(output.log).open("w", encoding="utf-8") as log_handle:
            subprocess.run(command, check=True, stdout=log_handle, stderr=log_handle)


rule aggregate_gene_trees:
    input:
        report=GENE_TREE_DIRECTORY_OUTPUTS["report"],
        treefile=GENE_TREE_DIRECTORY_OUTPUTS["treefile"],
        best_model_nex=GENE_TREE_DIRECTORY_OUTPUTS["best_model_nex"],
    output:
        manifest=GENE_TREE_MANIFEST,
        aggregate=GENE_TREE_AGGREGATE,
    params:
        support_mode=IQTREE_SUPPORT_MODE,
    shell:
        (
            "python3 -m scripts.finalize_gene_trees "
            "--best-model-nex {input.best_model_nex:q} "
            "--treefile {input.treefile:q} "
            "--report {input.report:q} "
            "--manifest {output.manifest:q} "
            "--output {output.aggregate:q} "
            "--support-mode {params.support_mode:q}"
        )
