import pandas as pd
from check_lipinski import process_lipinski
from fragments import process_fragments
from check_smiles_chemval import process_chemval

def process_pipeline(input_file, output_file):
    # Load input data
    input_df = pd.read_csv(input_file, names=["SMILES"])
    
    # Process Lipinski compliance
    lipinski_df = process_lipinski(input_df)
    
    # Add BRICS fragments
    fragments_df = process_fragments(lipinski_df)
    
    # Validate SMILES (chemval)
    final_df = process_chemval(fragments_df)
    
    # Save final results
    final_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process SMILES through Lipinski, fragments, and chemical validation.")
    parser.add_argument("--input", required=True, help="Input file containing SMILES strings.")
    parser.add_argument("--output", required=True, help="Output file to save results.")
    args = parser.parse_args()

    process_pipeline(args.input, args.output)
