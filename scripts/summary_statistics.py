import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def calculate_statistics(input_file, stats_output):
    # Load the dataset
    data = pd.read_csv(input_file)

    # Select numerical columns for analysis
    numerical_columns = ['HBA', 'HBD', 'LogP']

    # Calculate summary statistics
    summary = data[numerical_columns].describe()

    # Save statistics to a CSV file
    summary.to_csv(stats_output)

def generate_visualizations(input_file, plot_output):
    # Load the dataset
    data = pd.read_csv(input_file)

    # Select numerical columns for analysis
    numerical_columns = ['HBA', 'HBD', 'LogP']

    # Create a pairplot for numerical features
    sns.pairplot(data, vars=numerical_columns, diag_kind="kde", hue="Valid")

    # Save the visualization to a PNG file
    plt.savefig(plot_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate summary statistics and visualizations.")
    parser.add_argument('--input', required=True, help="Input CSV file")
    parser.add_argument('--stats-output', required=True, help="Output file for summary statistics")
    parser.add_argument('--plot-output', required=True, help="Output file for visualizations")

    args = parser.parse_args()

    calculate_statistics(args.input, args.stats_output)
    generate_visualizations(args.input, args.plot_output)
