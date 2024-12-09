What does this do?
- A [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to generate all possible isomers from a molecular formula ([Surge](https://github.com/StructureGenerator/surge))
- Checking the output (smiles strings) for the Lipinski rule of five ([rdkit.Lipinski](https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html))
- Generating synthons ([rdkit.BRICS](https://www.rdkit.org/docs/source/rdkit.Chem.BRICS.html))
- Some random analysis and vizualizations to populate the snakemake workflow.

How is this useful?
- It isn't

Requirements:
- Docker
- Snakemake

Usage:
- Make sure you have Docker and Snakemake installed
- Choose an input molecular formula in the surge_input.txt file
- Run the workflow with $ snakemake
- Outputs are stored in the output folder

Problems:
- There are no checks for the validity of molecules (e.g. geminal di- and triols, over constrained rings, unstable tautomers)