import time

# Generate a timestamp to append to outputs
timestamp = time.strftime("%Y%m%d_%H%M%S")

rule all:
    input:
        expand("lipinski_results_{timestamp}.csv", timestamp=timestamp)

rule build_pydev_image:
    output:
        "pydev_image.sentinel"
    shell:
        """
        if ! docker images | grep -q 'pydev'; then
            docker build -t pydev ./python-container
        fi
        touch {output}
        """

rule build_surge_image:
    output:
        "surge_image.sentinel"
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
        "surge_image.sentinel"
    output:
        "surge_output_{timestamp}.smi"
    params:
        timestamp=timestamp
    run:
        with open(input[0], 'r') as f:
            input_string = f.read().strip()
        shell(
            './surge.sh "{input_string}" {output}'
        )

rule check_lipinski:
    input:
        "surge_output_{timestamp}.smi",
        "pydev_image.sentinel"
    output:
        "lipinski_results_{timestamp}.csv"
    params:
        timestamp=timestamp
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python check_lipinski.py \
        --input {input[0]} --output {output}
        """
