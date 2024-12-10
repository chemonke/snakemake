# Load configuration
configfile: "config/config.yaml"

# Load dynamic paths
scripts = config["folders"]["scripts"]
output = config["folders"]["output"]
containers = config["folders"]["containers"]
sentinel = config["folders"]["sentinel"]


# Include modularized rules
include: "./rules/containers"
include: "./rules/data_generation"
include: "./rules/analytics"

rule all:
    input:
        expand(f"{output}/{{formula}}_surge_output.smi", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_lipinski_results_with_fragments.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_summary_statistics.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_pairplot.png", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_validation_error_stats.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_validation_error_visualization.png", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_fragment_analysis.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_fragment_frequency_plot.png", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_molecules.png", formula=lambda wildcards: get_formulas())


# Helper function to read formulas from formulas.txt
def get_formulas():
    try:
        with open("./config/surge_input.txt") as f:
            return [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        return []

    
