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
        # Outputs from the surge processing step
        expand(f"{output}/{{formula}}_surge_output.smi", formula=lambda wildcards: get_formulas()),
        # Outputs from Lipinski results processing
        expand(f"{output}/{{formula}}_lipinski_results_with_fragments.csv", formula=lambda wildcards: get_formulas()),
        # Outputs from statistical analysis
        expand(f"{output}/{{formula}}_summary_statistics.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_pairplot.png", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_validation_error_stats.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_validation_error_visualization.png", formula=lambda wildcards: get_formulas()),
        # Outputs from fragment analysis
        expand(f"{output}/{{formula}}_fragment_analysis.csv", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_fragment_frequency_plot.png", formula=lambda wildcards: get_formulas()),
        expand(f"{output}/{{formula}}_molecules.png", formula=lambda wildcards: get_formulas()),
        # Data insertion completion flags
        expand(f"{sentinel}/{{formula}}_data_inserted.flag", formula=lambda wildcards: get_formulas()),
        # Ensure Docker network and required containers are ready
        config["files"]["network_sentinel"],
        config["files"]["sql_sentinel"],
        config["files"]["pydev_sentinel"]




# Helper function to read formulas from formulas.txt
def get_formulas():
    try:
        with open("./config/surge_input.txt") as f:
            return [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        return []

    
