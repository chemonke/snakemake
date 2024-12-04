from rdkit import Chem
from rdkit.Chem import Draw, Lipinski, Descriptors
import os

def check_lipinski(smiles):
    """Check if a molecule complies with Lipinski's Rule of Five."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, True

    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)

    num_violations = sum([
        hbd > 5,
        hba > 10,
        mw > 500,
        logp > 5
    ])

    return num_violations <= 1, mol  # True if valid, False if violates

def generate_images(input_file, valid_output, invalid_output):
    """Generate images of valid and invalid molecules."""
    valid_mols = []
    invalid_mols = []

    # Ensure input file exists and is not empty
    if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
        raise FileNotFoundError(f"Input file '{input_file}' does not exist or is empty.")

    # Read SMILES strings and categorize molecules
    with open(input_file, 'r') as infile:
        for line in infile:
            smiles = line.strip()
            is_valid, mol = check_lipinski(smiles)
            if mol:
                if is_valid and len(valid_mols) < 5:
                    valid_mols.append(mol)
                elif not is_valid and len(invalid_mols) < 5:
                    invalid_mols.append(mol)
            if len(valid_mols) >= 5 and len(invalid_mols) >= 5:
                break

    # Save images to specified output files
    if valid_mols:
        Draw.MolsToImage(valid_mols).save(valid_output)
    if invalid_mols:
        Draw.MolsToImage(invalid_mols).save(invalid_output)

if __name__ == "__main__":
    import argparse

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Draw images of molecules based on Lipinski's rule compliance.")
    parser.add_argument("--input", required=True, help="Input file containing SMILES strings.")
    parser.add_argument("--valid-output", required=True, help="Output file for valid molecules image.")
    parser.add_argument("--invalid-output", required=True, help="Output file for invalid molecules image.")
    args = parser.parse_args()

    # Generate molecule images
    generate_images(args.input, args.valid_output, args.invalid_output)
