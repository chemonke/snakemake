rule all:
    input:
        "lipinski_results.csv"  # Final target of the workflow

rule build_pydev_image:
    output:
        "pydev_image.sentinel" # serves as a check for whether or not the image has been built
    shell:
        """
        if ! docker images | grep -q 'pydev'; then
            docker build -t pydev ./python-container
        fi
        touch {output}
        """

rule build_surge_image:
    output:
        "surge_image.sentinel" # serves as a check for whether or not the image has been built
    shell:
        """
        if ! docker images | grep -q 'surge'; then
            docker build -t surge ./surge
        fi
        touch {output}
        """

rule run_surge:
    input:
        "surge_input.txt",
        "surge_image.sentinel"  # Ensure the surge image is built
    output:
        "surge_output.smi"
    shell:
        """
        docker run --rm -v $(pwd):/workspace -w /workspace surge \
        ./surge.sh "$(cat {input[0]})" {output}
        """

rule check_lipinski:
    input:
        "surge_output.smi",
        "pydev_image.sentinel"  # Ensure the snakemake-dev image is built
    output:
        "lipinski_results.csv"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python check_lipinski.py \
        --input {input[0]} --output {output}
        """
