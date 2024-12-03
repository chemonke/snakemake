import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS
import argparse
import os

# Function to generate BRICS fragments
def get_brics_fragments(smiles):
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            fragments = list(BRICS.BRICSDecompose(molecule))
            return fragments  # Return the fragments as a list
        else:
            return ["Invalid SMILES"]  # Return as a list for consistency
    except Exception as e:
        return [f"Error: {e}"]  # Return as a list for consistency

# Process the CSV
def process_smiles_with_fragments(input_file, output_file):
    # Check if input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' does not exist.")

    df = pd.read_csv(input_file)
    
    # Generate BRICS fragments for each SMILES
    df['Fragments'] = df['SMILES'].apply(get_brics_fragments)
    
    # Save the updated DataFrame to a new file
    df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")


# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BRICS fragments for SMILES in a CSV file.")
    parser.add_argument("--input", required=True, help="Input CSV file containing SMILES strings.")
    parser.add_argument("--output", required=True, help="Output CSV file to save results.")
    args = parser.parse_args()

    # Run the processing function
    process_smiles_with_fragments(args.input, args.output)
