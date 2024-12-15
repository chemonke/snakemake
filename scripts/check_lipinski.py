from rdkit import Chem
from rdkit.Chem import Crippen, Lipinski, Descriptors
import pandas as pd

def check_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, True  # Return True for the Error column if invalid SMILES

    # Calculate Lipinski properties
    hbd = Lipinski.NumHDonors(mol)  # Hydrogen bond donors
    hba = Lipinski.NumHAcceptors(mol)  # Hydrogen bond acceptors
    mw = Descriptors.ExactMolWt(mol)  # Molecular weight
    logp = Crippen.MolLogP(mol)  # LogP

    # Check individual violations
    lipinski_violations = {
        "HBD > 5": hbd > 5,
        "HBA > 10": hba > 10,
        "MW > 500": mw > 500,
        "LogP > 5": logp > 5,
    }

    num_violations = sum(lipinski_violations.values())

    return {
        "MW": mw,
        "HBA": hba,
        "HBD": hbd,
        "LogP": logp,
        "lipinski_valid": num_violations <= 1,  # Valid if at most 1 violation
    }, False  # No error

def process_lipinski(df):
    results = []
    for smiles in df['SMILES']:
        properties, error = check_lipinski(smiles)
        if error:
            results.append({
                "SMILES": smiles,
                "MW": None,
                "HBA": None,
                "HBD": None,
                "LogP": None,
                "lipinski_valid": False,
                "Error": True,  # Set Error to True for invalid SMILES
            })
        else:
            results.append({
                "SMILES": smiles,
                **properties,
                "Error": False  # Set Error to False for valid SMILES
            })
    return pd.DataFrame(results)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Check Lipinski's Rule of Five compliance for SMILES.")
    parser.add_argument("--input", required=True, help="Input file containing SMILES strings.")
    parser.add_argument("--output", required=True, help="Output file to save results.")
    args = parser.parse_args()

    # Read input SMILES
    input_df = pd.read_csv(args.input, names=["SMILES"])
    output_df = process_lipinski(input_df)

    # Save results
    output_df.to_csv(args.output, index=False)
