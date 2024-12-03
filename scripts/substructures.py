import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import argparse
from collections import Counter

def analyze_fragments(input_file, output_file, plot_output):
    # Load the dataset
    print(f"Reading input file: {input_file}")
    data = pd.read_csv(input_file)
    
    # Flatten the list of fragments
    print("Processing fragments...")
    all_fragments = []
    for fragment_list in data['Fragments']:
        # Convert string representation of list to actual list
        fragments = eval(fragment_list) if isinstance(fragment_list, str) else []
        all_fragments.extend(fragments)
    
    # Count the occurrences of each fragment
    fragment_counts = Counter(all_fragments)
    
    # Convert to DataFrame and sort by frequency
    fragment_df = pd.DataFrame(fragment_counts.items(), columns=['Fragment', 'Count'])
    fragment_df = fragment_df.sort_values(by='Count', ascending=False)
    
    # Save the analysis results
    print(f"Saving fragment analysis to {output_file}")
    fragment_df.to_csv(output_file, index=False)
    
    # Draw the top 10 most common fragments
    print("Generating drawings for the most common fragments...")
    top_fragments = fragment_df.head(10)
    
    # Convert SMILES to RDKit molecule objects
    mols = [Chem.MolFromSmiles(smiles) for smiles in top_fragments['Fragment']]
    
    # Check for invalid molecules
    for idx, mol in enumerate(mols):
        if mol is None:
            print(f"Warning: Invalid SMILES string at index {idx}: {top_fragments.iloc[idx]['Fragment']}")
    
    # Draw molecules
    img = Draw.MolsToGridImage(
        mols, 
        legends=[f"{row['Fragment']} (Count: {row['Count']})" for _, row in top_fragments.iterrows()],
        molsPerRow=5, 
        subImgSize=(300, 300)
    )
    
    # Save the image
    print(f"Saving fragment drawings to {plot_output}")
    img.save(plot_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze and visualize fragment frequencies.")
    parser.add_argument('--input', required=True, help="Input CSV file")
    parser.add_argument('--output', required=True, help="Output CSV file for fragment analysis")
    parser.add_argument('--plot-output', required=True, help="Output PNG file for fragment frequency plot")

    args = parser.parse_args()
    
    analyze_fragments(args.input, args.output, args.plot_output)
