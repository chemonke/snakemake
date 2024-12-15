import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw
import argparse

def generate_images(csv_file, output_file, max_molecules=10):
    """Generate a single image showing valid and invalid molecules from precomputed CSV."""
    # Read the CSV file
    data = pd.read_csv(csv_file)

    # Ensure necessary columns are present
    required_columns = ["SMILES", "lipinski_valid", "chemval"]
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"CSV must contain '{col}' column.")

    # Filter data to only include molecules with chemval = True
    data = data[data["chemval"] == True]

    # Limit molecules to max_molecules per category
    valid_mols = [Chem.MolFromSmiles(smiles) for smiles in data[data["lipinski_valid"] == True]["SMILES"][:max_molecules] if Chem.MolFromSmiles(smiles)]
    invalid_mols = [Chem.MolFromSmiles(smiles) for smiles in data[data["lipinski_valid"] == False]["SMILES"][:max_molecules] if Chem.MolFromSmiles(smiles)]

    # Define grid settings
    sub_img_size = (200, 200)
    grid_size = (2, 5)  # 2 rows, 5 columns (adjust if needed)

    # Generate molecule images in grids
    valid_img = Draw.MolsToGridImage(valid_mols, molsPerRow=grid_size[1], subImgSize=sub_img_size) if valid_mols else None
    invalid_img = Draw.MolsToGridImage(invalid_mols, molsPerRow=grid_size[1], subImgSize=sub_img_size) if invalid_mols else None

    # Combine into a single image
    labels = ["Valid Molecules", "Invalid Molecules"]
    images = [valid_img, invalid_img]
    combined_height = (sub_img_size[1] * grid_size[0] + 50) * len(images)  # Include space for labels
    combined_img = Image.new("RGB", (sub_img_size[0] * grid_size[1], combined_height), "white")
    draw = ImageDraw.Draw(combined_img)
    font = ImageFont.load_default()

    # Paste each grid into the combined image
    y_offset = 0
    for img, label in zip(images, labels):
        if img:
            combined_img.paste(img, (0, y_offset))
        draw.text((10, y_offset + sub_img_size[1] * grid_size[0] + 10), label, fill="black", font=font)
        y_offset += sub_img_size[1] * grid_size[0] + 50  # Move down for the next grid + label

    # Save the combined image
    combined_img.save(output_file)


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Draw images of molecules based on precomputed Lipinski results.")
    parser.add_argument("--csv", required=True, help="Input CSV file containing SMILES strings and Lipinski compliance.")
    parser.add_argument("--output", required=True, help="Output file for the combined molecules image.")
    parser.add_argument("--max-molecules", type=int, default=10, help="Maximum number of molecules to display per category.")
    args = parser.parse_args()

    # Generate molecule images
    generate_images(args.csv, args.output, max_molecules=args.max_molecules)
