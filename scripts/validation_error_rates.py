import pandas as pd
import matplotlib.pyplot as plt

def visualize_data(input_file, stats_output, plot_output):
    # Load the file
    data = pd.read_csv(input_file)

    # Ensure necessary columns are present
    required_columns = ['chemval', 'lipinski_valid']
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"Input CSV must contain '{col}' column.")

    # Categorize data into three groups
    chem_invalid_count = len(data[data['chemval'] == False])
    chem_and_lipinski_valid_count = len(data[(data['chemval'] == True) & (data['lipinski_valid'] == True)])
    chem_valid_lipinski_invalid_count = len(data[(data['chemval'] == True) & (data['lipinski_valid'] == False)])

    # Calculate total and validation percentage
    total_count = len(data)
    validation_percentage = round((chem_and_lipinski_valid_count / total_count) * 100,2) if total_count > 0 else 0

    # Prepare statistics for output
    stats_df = pd.DataFrame([{
        "chemically_invalid": chem_invalid_count,
        "chemically_valid_and_lipinski_valid": chem_and_lipinski_valid_count,
        "chemically_valid_but_lipinski_invalid": chem_valid_lipinski_invalid_count,
        "total": total_count,
        "validation_percentage": validation_percentage
    }])
    stats_df.to_csv(stats_output, index=False)

    # Prepare data for the bar chart
    bar_data = pd.DataFrame({
        "Type": [
            "Chemically Invalid",
            "Chemically valid and Lipinski Valid",
            "Chemically valid but Lipinski Invalid"
        ],
        "Count": [
            chem_invalid_count,
            chem_and_lipinski_valid_count,
            chem_valid_lipinski_invalid_count
        ]
    })

    # Create the bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(bar_data['Type'], bar_data['Count'], color=['lightcoral', 'skyblue', 'gold'], edgecolor='k', alpha=0.7)
    
    # Add counts above bars for better readability
    for i, value in enumerate(bar_data['Count']):
        plt.text(i, value + 5, str(value), ha='center', fontsize=10)

    plt.title('Distribution of Chemical and Lipinski Validity', fontsize=14)
    plt.xlabel('Validation Category', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xticks(rotation=15, fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(plot_output)
    plt.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate bar chart for chemical and Lipinski validation counts.")
    parser.add_argument('--input', required=True, help="Input CSV file")
    parser.add_argument('--stats-output', required=True, help="Output file for statistics")
    parser.add_argument('--plot-output', required=True, help="Output file for bar chart")
    args = parser.parse_args()

    visualize_data(args.input, args.stats_output, args.plot_output)
