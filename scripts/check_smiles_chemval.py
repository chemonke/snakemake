import argparse
import csv
from rdkit import Chem
from rdkit.Chem import FilterCatalog

# Initialize PAINS filter catalog
def initialize_pains_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog.FilterCatalog(params)

catalog = initialize_pains_catalog()

# Define validation function
def check_violations(smiles):
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

    # chemval is True if no violations, False otherwise
    chemval = len(violations) == 0

    return chemval, violations

# CSV Processing Function
def process_csv(input_file, output_file):
    with open(input_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ["chemval", "chem_violation"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            smiles = row.get("SMILES", "")
            chemval, violations = check_violations(smiles)
            row["chemval"] = chemval
            row["chem_violation"] = ", ".join(violations) if violations else ""
            writer.writerow(row)

# Main Entry Point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate SMILES in a CSV file using PAINS filters and additional rules. Add chemval and chem_violation columns.")
    parser.add_argument("--input", required=True, help="Path to input CSV file.")
    parser.add_argument("--output", required=True, help="Path to output CSV file.")
    args = parser.parse_args()

    process_csv(args.input, args.output)
