import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS

def get_brics_fragments(smiles):
    """Generate BRICS fragments for a given SMILES."""
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            return list(BRICS.BRICSDecompose(molecule))  # Return fragments as a list
        else:
            return ["Invalid SMILES"]  # Return as a list for consistency
    except Exception as e:
        return [f"Error: {e}"]  # Return as a list for consistency

def process_fragments(df):
    """
    Generate BRICS fragments for SMILES strings in a DataFrame.
    
    Args:
        df (pd.DataFrame): Input DataFrame with a 'SMILES' column.
        
    Returns:
        pd.DataFrame: Updated DataFrame with a 'Fragments' column.
    """
    df = df.copy()  # Avoid modifying the input DataFrame
    df['Fragments'] = df['SMILES'].apply(get_brics_fragments)
    return df

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Generate BRICS fragments for SMILES in a CSV file.")
    parser.add_argument("--input", required=True, help="Input CSV file containing SMILES strings.")
    parser.add_argument("--output", required=True, help="Output CSV file to save results.")
    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file '{args.input}' does not exist.")

    # Load the input DataFrame
    input_df = pd.read_csv(args.input)

    # Generate fragments
    output_df = process_fragments(input_df)

    # Save the updated DataFrame
    output_df.to_csv(args.output, index=False)
    print(f"Processed data saved to {args.output}")
