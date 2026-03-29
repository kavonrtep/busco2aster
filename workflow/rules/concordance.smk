rule infer_gene_concordance:
    input:
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        gene_trees=GENE_TREE_AGGREGATE,
    output:
        command=GCF_OUTPUTS["command"],
        stat=GCF_OUTPUTS["stat"],
        tree=GCF_OUTPUTS["tree"],
        branch=GCF_OUTPUTS["branch"],
        log=GCF_OUTPUTS["log"],
        completion=GCF_OUTPUTS["completion"],
    params:
        executable=config["iqtree_executable"],
        prefix=GCF_OUTPUTS["prefix"],
    threads:
        get_thread_count("concordance")
    run:
        import subprocess
        from pathlib import Path

        from scripts.concordance import (
            build_iqtree_gcf_command,
            write_concordance_command_script,
        )

        command = build_iqtree_gcf_command(
            executable=params.executable,
            reference_tree_path=input.species_tree,
            gene_tree_path=input.gene_trees,
            prefix=params.prefix,
            threads=threads,
        )
        write_concordance_command_script(Path(output.command), command)
        subprocess.run(command, check=True)
        Path(output.completion).touch()
