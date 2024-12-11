rule run_surge:
    input:
        config["files"]["surge_sentinel"],
        config["files"]["pydev_sentinel"]
    output:
        f"{output}/{{formula}}_surge_output.smi"
    params:
        formula=lambda wildcards: wildcards.formula
    shell:
        """
        ./scripts/surge.sh "{params.formula}" {output}
        """