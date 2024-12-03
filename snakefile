rule all:
    input:
        "lipinski_results.csv"  # Final target of the workflow

rule run_surge:
    input:
        "surge_input.txt"
    output:
        "surge_output.smi"
    run:
        with open(input[0], 'r') as f:
            input_string = f.read().strip()  # Reading the string from the txt file
        shell(
            './surge.sh "{input_string}" {output[0]}'
        )

rule check_lipinski:
    input:
        "surge_output.smi"
    output:
        "lipinski_results.csv"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace snakemake-dev \
        /opt/conda/envs/snakemake-env/bin/python check_lipinski.py \
        --input {input} --output {output}
        """
