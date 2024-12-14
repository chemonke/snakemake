import pandas as pd
import matplotlib.pyplot as plt

def visualize_data(input_file, stats_output, plot_output):
    # Load the file
    data = pd.read_csv(input_file)

    # Count valid and invalid entries
    valid_count = data['Valid'].apply(lambda x: str(x).lower() == 'true').sum() if 'Valid' in data.columns else 0
    invalid_count = len(data) - valid_count

    # Calculate the validation percentage
    total_count = valid_count + invalid_count
    validation_percentage = (valid_count / total_count) * 100 if total_count > 0 else 0

    # Create a DataFrame with the counts and percentage
    stats_df = pd.DataFrame([{
        "valid": valid_count,
        "invalid": invalid_count,
        "validation_percentage": validation_percentage
    }])
    stats_df.to_csv(stats_output, index=False)

    # Prepare data for the bar chart (only valid and invalid counts)
    bar_data = pd.DataFrame({
        "Type": ["valid", "invalid"],
        "Count": [valid_count, invalid_count]
    })

    # Create the bar chart
    plt.figure(figsize=(8, 6))
    plt.bar(bar_data['Type'], bar_data['Count'], color=['skyblue', 'salmon'], edgecolor='k', alpha=0.7)
    
    # Add counts above bars for better readability
    for i, value in enumerate(bar_data['Count']):
        plt.text(i, value + 5, str(value), ha='center', fontsize=10)

    plt.title('Distribution of Outputs', fontsize=14)
    plt.xlabel('Type', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xticks(rotation=0)
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
