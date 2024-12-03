import pandas as pd
import matplotlib.pyplot as plt
import argparse

def calculate_rates(input_file, stats_output):
    # Load the dataset
    data = pd.read_csv(input_file)
    
    # Calculate proportions
    valid_counts = data['Valid'].value_counts(normalize=True) * 100
    error_counts = data['Error'].value_counts(normalize=True) * 100
    
    # Combine results
    results = pd.DataFrame({
        'Category': ['Valid True', 'Valid False', 'Error False', 'Error True'],
        'Proportion (%)': [
            valid_counts.get(True, 0),
            valid_counts.get(False, 0),
            error_counts.get(False, 0),
            error_counts.get(True, 0)
        ]
    })
    
    # Save results
    results.to_csv(stats_output, index=False)

def visualize_rates(input_file, plot_output):
    print(f"Generating visualizations from {input_file}")
    data = pd.read_csv(input_file)
    
    valid_counts = data['Valid'].value_counts(normalize=True)
    error_counts = data['Error'].value_counts(normalize=True)
    
    # Ensure labels match the actual data
    valid_labels = [f"{label}" for label in valid_counts.index]
    error_labels = [f"{label}" for label in error_counts.index]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Valid Plot
    axes[0].pie(valid_counts, labels=valid_labels, autopct='%1.1f%%', colors=['#4CAF50', '#F44336'][:len(valid_counts)])
    axes[0].set_title('Proportion of Valid Molecules')
    
    # Error Plot
    axes[1].pie(error_counts, labels=error_labels, autopct='%1.1f%%', colors=['#4CAF50', '#F44336'][:len(error_counts)])
    axes[1].set_title('Proportion of Errors')
    
    print(f"Saving plots to {plot_output}")
    plt.tight_layout()
    plt.savefig(plot_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute validation and error rates.")
    parser.add_argument('--input', required=True, help="Input CSV file")
    parser.add_argument('--stats-output', required=True, help="Output file for statistics")
    parser.add_argument('--plot-output', required=True, help="Output file for visualization")

    args = parser.parse_args()

    calculate_rates(args.input, args.stats_output)
    visualize_rates(args.input, args.plot_output)
