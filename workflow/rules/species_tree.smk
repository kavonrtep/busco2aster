rule prepare_wastral_gene_trees:
    input:
        gene_trees=GENE_TREE_AGGREGATE,
    output:
        WASTRAL_GENE_TREE_INPUT,
    params:
        support_mode=IQTREE_SUPPORT_MODE,
    run:
        from pathlib import Path

        from scripts.species_tree import prepare_wastral_gene_tree_input

        prepare_wastral_gene_tree_input(
            input_path=Path(input.gene_trees),
            output_path=Path(output[0]),
            support_mode=params.support_mode,
        )


rule infer_species_tree_wastral:
    input:
        gene_trees=WASTRAL_GENE_TREE_INPUT,
    output:
        command=WASTRAL_OUTPUTS["command"],
        treefile=WASTRAL_OUTPUTS["treefile"],
        log=WASTRAL_OUTPUTS["log"],
        completion=WASTRAL_OUTPUTS["completion"],
    params:
        executable=WASTRAL_EXECUTABLE,
        support_mode=IQTREE_SUPPORT_MODE,
    threads:
        get_thread_count("species_tree")
    run:
        import subprocess
        from pathlib import Path

        from scripts.species_tree import build_aster_command, write_species_tree_command_script

        command = build_aster_command(
            backend="wastral",
            executable=params.executable,
            input_path=input.gene_trees,
            output_path=output.treefile,
            threads=threads,
            support_mode=params.support_mode,
        )
        write_species_tree_command_script(Path(output.command), command, Path(output.log))
        with Path(output.log).open("w", encoding="utf-8") as log_handle:
            subprocess.run(command, check=True, stderr=log_handle)
        Path(output.completion).touch()


rule infer_species_tree_astral4:
    input:
        gene_trees=GENE_TREE_AGGREGATE,
    output:
        command=ASTRAL4_OUTPUTS["command"],
        treefile=ASTRAL4_OUTPUTS["treefile"],
        log=ASTRAL4_OUTPUTS["log"],
        completion=ASTRAL4_OUTPUTS["completion"],
    params:
        executable=ASTRAL4_EXECUTABLE,
    threads:
        get_thread_count("species_tree")
    run:
        import subprocess
        from pathlib import Path

        from scripts.species_tree import build_aster_command, write_species_tree_command_script

        command = build_aster_command(
            backend="astral4",
            executable=params.executable,
            input_path=input.gene_trees,
            output_path=output.treefile,
            threads=threads,
        )
        write_species_tree_command_script(Path(output.command), command, Path(output.log))
        with Path(output.log).open("w", encoding="utf-8") as log_handle:
            subprocess.run(command, check=True, stderr=log_handle)
        Path(output.completion).touch()
