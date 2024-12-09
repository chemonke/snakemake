import pandas as pd
import matplotlib.pyplot as plt

def visualize_data(input_file, stats_output, plot_output):
    # Load the file
    data = pd.read_csv(input_file)

    # Count True/False in Valid column
    if 'Valid' in data.columns:
        valid_counts = data['Valid'].value_counts()
        valid_counts.to_csv(stats_output, header=['Count'])

        # Create the bar chart for Valid column
        plt.figure(figsize=(8, 6))
        valid_counts.plot(kind='bar', color=['skyblue', 'salmon'], edgecolor='k', alpha=0.7)
        plt.title('Distribution of Outputs', fontsize=14)
        plt.xlabel('Valid', fontsize=12)
        plt.ylabel('Count', fontsize=12)
        plt.xticks(ticks=[0, 1], labels=['True', 'False'], rotation=0)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.savefig(plot_output)
        plt.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate bar chart for Valid counts.")
    parser.add_argument('--input', required=True, help="Input CSV file")
    parser.add_argument('--stats-output', required=True, help="Output file for statistics")
    parser.add_argument('--plot-output', required=True, help="Output file for bar chart")
    args = parser.parse_args()

    visualize_data(args.input, args.stats_output, args.plot_output)
