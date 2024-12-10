rule check_lipinski:
    input:
        f"{output}/{{formula}}_surge_output.smi",
        config["files"]["pydev_sentinel"]
    output:
        f"{output}/{{formula}}_lipinski_results.csv"
    shell:
        """
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
/opt/conda/envs/pydev/bin/python {scripts}/check_lipinski.py \
--input {input[0]} --output {output}
        """

rule add_fragments:
    input:
        f"{output}/{{formula}}_lipinski_results.csv"
    output:
        f"{output}/{{formula}}_lipinski_results_with_fragments.csv"
    shell:
        """
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
/opt/conda/envs/pydev/bin/python {scripts}/fragments.py \
--input {input} --output {output}
        """

rule generate_summary_statistics:
    input:
        f"{output}/{{formula}}_lipinski_results_with_fragments.csv"
    output:
        f"{output}/{{formula}}_summary_statistics.csv",
        f"{output}/{{formula}}_pairplot.png"
    shell:
        """
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
/opt/conda/envs/pydev/bin/python {scripts}/summary_statistics.py \
--input {input} \
--stats-output {output[0]} \
--plot-output {output[1]}
        """

rule draw_molecules:
    input:
        f"{output}/{{formula}}_lipinski_results.csv"
    output:
        f"{output}/{{formula}}_molecules.png"
    shell:
        """
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
/opt/conda/envs/pydev/bin/python {scripts}/draw_molecules.py \
--csv {input} --output {output} --max-molecules 10
        """

rule fragment_analysis:
    input:
        f"{output}/{{formula}}_lipinski_results_with_fragments.csv"
    output:
        f"{output}/{{formula}}_fragment_analysis.csv",
        f"{output}/{{formula}}_fragment_frequency_plot.png"
    shell:
        """
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
/opt/conda/envs/pydev/bin/python {scripts}/substructures.py \
--input {input} --output {output[0]} --plot-output {output[1]}
        """

rule validation_error_analysis:
    input:
        f"{output}/{{formula}}_lipinski_results_with_fragments.csv"
    output:
        f"{output}/{{formula}}_validation_error_stats.csv",
        f"{output}/{{formula}}_validation_error_visualization.png"
    shell:
        """
docker run --rm --user $(id -u):$(id -g) -v $(pwd):/workspace -w /workspace pydev \
/opt/conda/envs/pydev/bin/python {scripts}/validation_error_rates.py \
--input {input} \
--stats-output {output[0]} \
--plot-output {output[1]}
        """