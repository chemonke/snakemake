import yaml

# Load config.yaml
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Use folder structure from config
surge_input = config["files"]["surge_input"]

# Process the surge_input file
try:
    with open(surge_input, "r") as f:
        formulas = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"Error: File {surge_input} not found.")
    formulas = []

# Output the formulas
for formula in formulas:
    print(formula)
