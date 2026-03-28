rule summarize_busco:
    input:
        manifest=ancient(VALIDATED_MANIFEST),
        busco_artifacts=BUSCO_SUMMARY_INPUTS,
    output:
        summary=BUSCO_SUMMARY_TABLE,
        records=BUSCO_RECORDS_TABLE,
    shell:
        (
            "python3 -m scripts.summarize_busco "
            "--manifest {input.manifest:q} "
            "--summary-output {output.summary:q} "
            "--records-output {output.records:q}"
        )
