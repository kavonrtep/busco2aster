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
        executable=IQTREE_EXECUTABLE,
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


rule infer_site_concordance:
    input:
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        alignments_complete=ALIGNMENTS_COMPLETE,
    output:
        command=SCFL_OUTPUTS["command"],
        stat=SCFL_OUTPUTS["stat"],
        tree=SCFL_OUTPUTS["tree"],
        branch=SCFL_OUTPUTS["branch"],
        log=SCFL_OUTPUTS["log"],
        completion=SCFL_OUTPUTS["completion"],
    params:
        executable=IQTREE_EXECUTABLE,
        prefix=SCFL_OUTPUTS["prefix"],
        alignment_dir=ALIGNMENT_DIR,
        quartets=IQTREE_SCFL_QUARTETS,
        seqtype="AA",
        model=IQTREE_SCFL_MODEL,
    threads:
        get_thread_count("concordance")
    run:
        import subprocess
        from pathlib import Path

        from scripts.concordance import (
            build_iqtree_scfl_command,
            write_concordance_command_script,
        )

        command = build_iqtree_scfl_command(
            executable=params.executable,
            reference_tree_path=input.species_tree,
            alignment_dir=params.alignment_dir,
            prefix=params.prefix,
            threads=threads,
            quartets=params.quartets,
            seqtype=params.seqtype,
            model=params.model,
        )
        write_concordance_command_script(Path(output.command), command)
        subprocess.run(command, check=True)
        Path(output.completion).touch()


rule annotate_species_tree_quartets:
    input:
        species_tree=DEFAULT_SPECIES_TREE_OUTPUTS["treefile"],
        gene_trees=WASTRAL_GENE_TREE_INPUT,
    output:
        command=WASTRAL_QUARTET_OUTPUTS["command"],
        tree=WASTRAL_QUARTET_OUTPUTS["tree"],
        freqquad=WASTRAL_QUARTET_OUTPUTS["freqquad"],
        log=WASTRAL_QUARTET_OUTPUTS["log"],
        completion=WASTRAL_QUARTET_OUTPUTS["completion"],
    params:
        executable=WASTRAL_EXECUTABLE,
        support_mode=IQTREE_SUPPORT_MODE,
        annotation_mode=3,
    threads:
        get_thread_count("concordance")
    run:
        import shutil
        import subprocess
        import tempfile
        from pathlib import Path

        from scripts.species_tree import build_aster_command, write_species_tree_command_script

        work_root = Path("work") / "concordance"
        work_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="wastral_quartets.", dir=work_root.as_posix()) as tmpdir:
            tmp_path = Path(tmpdir).resolve()
            tree_output = tmp_path / "annotated.tre"
            freqquad_output = tmp_path / "freqQuad.csv"
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
            if not freqquad_output.is_file():
                raise FileNotFoundError(f"wASTRAL quartet scoring did not write expected file: {freqquad_output}")
            Path(output.tree).parent.mkdir(parents=True, exist_ok=True)
            shutil.move(tree_output.as_posix(), output.tree)
            shutil.move(freqquad_output.as_posix(), output.freqquad)
        Path(output.completion).touch()
