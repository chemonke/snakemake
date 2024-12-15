# Load configuration
configfile: "config/config.yaml"

# Load dynamic paths
scripts = config["folders"]["scripts"]
output = config["folders"]["output"]
containers = config["folders"]["containers"]
sentinel = config["folders"]["sentinel"]


# Include modularized rules
include: "./rules/containers.smk"
include: "./rules/data_generation.smk"
include: "./rules/analytics.smk"
include: "./rules/sql.smk"

rule all:
    input:
        # Ensure Docker containers and network are ready
        config["files"]["network_sentinel"],
        config["files"]["sql_sentinel"],
        config["files"]["pydev_sentinel"],

        # Surge processing outputs
        expand(f"{output}/{{formula}}_surge_output.smi", formula=lambda wildcards: get_formulas()),

        # Lipinski results processing outputs
        expand(f"{output}/{{formula}}_results.csv", formula=lambda wildcards: get_formulas()),

        # Outputs from statistical analysis
        expand(f"{output}/{{formula}}_summary_statistics.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_pairplot.png", formula=lambda wildcards: get_formulas()),

        # Validation error analysis
        expand(f"{output}/{{formula}}_validation_error_stats.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_validation_error_visualization.png", formula=lambda wildcards: get_formulas()),

        # Fragment analysis outputs
        expand(f"{output}/{{formula}}_fragment_analysis.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_fragment_frequency_plot.png", formula=lambda wildcards: get_formulas()),

        # Molecule visualization
        expand(f"{output}/{{formula}}_molecules.png", formula=lambda wildcards: get_formulas()),

        # SQL insertion outputs
        expand(f"{sentinel}/{{formula}}_data_inserted.flag", formula=lambda wildcards: get_formulas())





# Helper function to read formulas from formulas.txt
def get_formulas():
    try:
        with open("./config/surge_input.txt") as f:
            return [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        return []