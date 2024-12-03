import time

# Generate a timestamp to append to outputs
timestamp = time.strftime("%Y%m%d_%H%M%S")

rule all:
    input:
        f"output/lipinski_results_with_fragments_{timestamp}.csv",
        f"output/summary_statistics_{timestamp}.csv",
        f"output/visualizations_{timestamp}.png"


rule build_pydev_image:
    output:
        "sentinel/pydev_image.sentinel"
    shell:
        """
        if ! docker images | grep -q 'pydev'; then
            docker build -t pydev ./python-container
        fi
        touch {output}
        """

rule build_surge_image:
    output:
        "sentinel/surge_image.sentinel"
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
        "sentinel/surge_image.sentinel"
    output:
        "output/surge_output_{timestamp}.smi"
    params:
        timestamp=timestamp
    run:
        with open(input[0], 'r') as f:
            input_string = f.read().strip()
        shell(
            './scripts/surge.sh "{input_string}" {output}'
        )


rule check_lipinski:
    input:
        "output/surge_output_{timestamp}.smi",
        "sentinel/pydev_image.sentinel"
    output:
        "output/lipinski_results_{timestamp}.csv"
    params:
        timestamp=timestamp
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/check_lipinski.py \
        --input {input[0]} --output {output}
        """

rule add_fragments:
    input:
        f"output/lipinski_results_{timestamp}.csv"
    output:
        f"output/lipinski_results_with_fragments_{timestamp}.csv"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/fragments.py \
        --input {input} --output {output}
        """

rule generate_summary_statistics:
    input:
        f"output/lipinski_results_with_fragments_{timestamp}.csv"
    output:
        f"output/summary_statistics_{timestamp}.csv",
        f"output/visualizations_{timestamp}.png"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/summary_statistics.py \
        --input {input} --stats-output {output[0]} --plot-output {output[1]}
        """
