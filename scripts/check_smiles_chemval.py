import os
import csv
from rdkit import Chem
from rdkit.Chem import FilterCatalog

def initialize_pains_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog.FilterCatalog(params)

catalog = initialize_pains_catalog()

def check_violations(smiles):
    """Check chemical rules and return validation status and violations."""
    violations = []

    # Convert SMILES to RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, ["Invalid SMILES"]

    # Check for PAINS violations
    if catalog.GetFirstMatch(mol):
        violations.append("PAINS")

    # Rule: No triple bonds in rings n<7
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) < 7:
            submol = Chem.PathToSubmol(mol, ring)
            if submol.HasSubstructMatch(Chem.MolFromSmarts("C#C")):
                violations.append("Triple bond in small ring")

    # Rule: No allenes
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C=C")):
        violations.append("Allene")

    # Rule: No geminal diols/triols
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[C](O)(O)")) or mol.HasSubstructMatch(Chem.MolFromSmarts("[C](O)(O)(O)")):
        violations.append("Geminal diol/triol")

    # Rule: No enols
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C]-[O]")):
        violations.append("Enol")

    # Rule: No azides
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[N]=[N+]=[N-]")):
        violations.append("Azide")

    # Rule: No hydrazines
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[N][N]")):
        violations.append("Hydrazine")

    # Rule: No peroxides
    if mol.HasSubstructMatch(Chem.MolFromSmarts("OO")):
        violations.append("Peroxide")

    # chemval is True if no violations, False otherwise
    chemval = len(violations) == 0

    return chemval, violations

def process_chemval(df):
    """Process a DataFrame to validate SMILES and append chemval columns."""
    df = df.copy()
    chemval_results = df['SMILES'].apply(check_violations)
    df['chemval'] = chemval_results.apply(lambda x: x[0])
    
    # Replace NaN or empty values in 'chem_violation' with "N/A"
    df['chem_violation'] = chemval_results.apply(lambda x: ", ".join(x[1]) if x[1] else "N/A")
    
    return df

if __name__ == "__main__":
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser(description="Validate SMILES in a CSV file using PAINS filters and additional rules.")
    parser.add_argument("--input", required=True, help="Path to input CSV file.")
    parser.add_argument("--output", required=True, help="Path to output CSV file.")
    args = parser.parse_args()

    # Load input CSV
    input_df = pd.read_csv(args.input)

    # Process chemical validation
    output_df = process_chemval(input_df)

    # Save results
    output_df.to_csv(args.output, index=False)
    print(f"Processed data saved to {args.output}")
