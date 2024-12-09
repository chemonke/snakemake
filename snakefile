formulas = shell("""
    docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
    /opt/conda/envs/pydev/bin/python ./scripts/read_formulae.py
""", read=True).strip().split("\n")

rule all:
    input:
        expand("{formula}_lipinski_results_with_fragments.csv", formula=formulas),
        expand("{formula}_summary_statistics.csv", formula=formulas),
        expand("{formula}_pairplot.png", formula=formulas),
        expand("{formula}_validation_error_stats.csv", formula=formulas),
        expand("{formula}_validation_error_visualization.png", formula=formulas),
        expand("{formula}_fragment_analysis.csv", formula=formulas),
        expand("{formula}_fragment_frequency_plot.png", formula=formulas),
        expand("{formula}_molecules.png", formula=formulas)

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
        "sentinel/surge_image.sentinel"
    output:
        "{formula}_surge_output.smi"
    params:
        formula="{formula}"
    shell:
        """
        ./scripts/surge.sh "{params.formula}" {output}
        """

rule check_lipinski:
    input:
        "{formula}_surge_output.smi",
        "sentinel/pydev_image.sentinel"
    output:
        "{formula}_lipinski_results.csv"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/check_lipinski.py \
        --input {input[0]} --output {output}
        """

rule add_fragments:
    input:
        "{formula}_lipinski_results.csv"
    output:
        "{formula}_lipinski_results_with_fragments.csv"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/fragments.py \
        --input {input} --output {output}
        """

rule generate_summary_statistics:
    input:
        "{formula}_lipinski_results_with_fragments.csv"
    output:
        "{formula}_summary_statistics.csv",
        "{formula}_pairplot.png"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/summary_statistics.py \
        --input {input} --stats-output {output[0]} --plot-output {output[1]}
        """

rule validate_error_rates:
    input:
        "{formula}_lipinski_results_with_fragments.csv"
    output:
        "{formula}_validation_error_stats.csv",
        "{formula}_validation_error_visualization.png"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace \
        -e MPLCONFIGDIR=/workspace/.config/matplotlib pydev \
        /opt/conda/envs/pydev/bin/python scripts/validation_error_rates.py \
        --input {input} --stats-output {output[0]} --plot-output {output[1]}
        """

rule analyze_fragments:
    input:
        "{formula}_lipinski_results_with_fragments.csv"
    output:
        "{formula}_fragment_analysis.csv",
        "{formula}_fragment_frequency_plot.png"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace \
        -e MPLCONFIGDIR=/workspace/.config/matplotlib pydev \
        /opt/conda/envs/pydev/bin/python scripts/substructures.py \
        --input {input} --output {output[0]} --plot-output {output[1]}
        """

rule draw_molecules:
    input:
        "{formula}_lipinski_results.csv"
    output:
        "{formula}_molecules.png"
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/draw_molecules.py \
        --csv {input} --output {output} --max-molecules 10
        """
