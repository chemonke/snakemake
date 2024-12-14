from rdkit import Chem
from rdkit.Chem import Crippen, Lipinski, Descriptors
import csv

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


def process_surge_output(input_file, output_file):
    results = []
    with open(input_file, 'r') as infile:
        for line in infile:
            smiles = line.strip()
            properties, error = check_lipinski(smiles)
            if error:
                results.append({
                    "SMILES": smiles,
                    "MW": None,
                    "HBA": None,
                    "HBD": None,
                    "LogP": None,
                    "Valid": False,
                    "Error": True,  # Set Error to True for invalid SMILES
                })
            else:
                results.append({
                    "SMILES": smiles,
                    **properties,
                    "Error": False  # Set Error to False for valid SMILES
                })

    # Explicitly set the desired column order
    fieldnames = ["SMILES", "MW", "HBA", "HBD", "LogP", "lipinski_valid", "Error"]

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Check Lipinski's Rule of Five compliance for SMILES.")
    parser.add_argument("--input", required=True, help="Input file containing SMILES strings.")
    parser.add_argument("--output", required=True, help="Output file to save results.")
    args = parser.parse_args()

    process_surge_output(args.input, args.output)
