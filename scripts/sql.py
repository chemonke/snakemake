import pymysql
import pandas as pd
import json  # For JSON validation
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Insert data into MariaDB")
parser.add_argument("--input", required=True, help="Path to the input CSV file")
parser.add_argument("--formula", required=True, help="Chemical formula being processed")
parser.add_argument("--user", required=True, help="Database username")
parser.add_argument("--password", required=True, help="Database password")
parser.add_argument("--host", default="mariadb", help="Database host")
parser.add_argument("--database", required=True, help="Database name")
args = parser.parse_args()

# Connect to the database
try:
    conn = pymysql.connect(
        host=args.host,
        user=args.user,
        password=args.password,
        database=args.database
    )
    cursor = conn.cursor()
    print("Successfully connected to the database.")
except pymysql.MySQLError as e:
    print(f"Error connecting to the database: {e}")
    exit(1)

# Create parent table for chemical formulae (without `description`)
create_formulae_table_query = """
CREATE TABLE IF NOT EXISTS chemical_formulae (
    id INT AUTO_INCREMENT PRIMARY KEY,
    formula TEXT NOT NULL,
    mw FLOAT DEFAULT NULL,
    lipinski_valid_percent FLOAT DEFAULT NULL,
    UNIQUE KEY (formula(255))
);
"""

# Create child table for lipinski_results
create_results_table_query = """
CREATE TABLE IF NOT EXISTS lipinski_results (
    id INT AUTO_INCREMENT PRIMARY KEY,
    formula_id INT NOT NULL,
    smiles TEXT NOT NULL,
    hba INT,
    hbd INT,
    logp FLOAT,
    lipinski_valid BOOLEAN,
    error BOOLEAN,
    fragments JSON,
    FOREIGN KEY (formula_id) REFERENCES chemical_formulae(id) ON DELETE CASCADE,
    UNIQUE KEY (formula_id, smiles(255))
);
"""

# Execute table creation queries
try:
    cursor.execute(create_formulae_table_query)
    print("Table `chemical_formulae` is ready.")
    cursor.execute(create_results_table_query)
    print("Table `lipinski_results` is ready.")
except pymysql.MySQLError as e:
    print(f"Error creating tables: {e}")
    cursor.close()
    conn.close()
    exit(1)

# Load the main CSV into a DataFrame
try:
    df = pd.read_csv(args.input)
    print(f"Successfully loaded data from {args.input}.")
except FileNotFoundError:
    print(f"File not found: {args.input}")
    cursor.close()
    conn.close()
    exit(1)
except pd.errors.ParserError as e:
    print(f"Error parsing CSV file: {e}")
    cursor.close()
    conn.close()
    exit(1)

# Calculate molecular weight (MW) as the average MW from the CSV
computed_mw = df["MW"].mean() if "MW" in df.columns else None

# Insert the formula into `chemical_formulae` with MW
formula_query = """
INSERT INTO chemical_formulae (formula, mw)
VALUES (%s, %s)
ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), mw=VALUES(mw);
"""
try:
    cursor.execute(formula_query, (args.formula, computed_mw))
    formula_id = cursor.lastrowid
    print(f"Chemical formula `{args.formula}` recorded with ID: {formula_id} and MW: {computed_mw}")
except pymysql.MySQLError as e:
    print(f"Error recording chemical formula: {e}")
    cursor.close()
    conn.close()
    exit(1)

# Load the validation percentage from the stats CSV
stats_file = f"output/{args.formula}_validation_error_stats.csv"
try:
    stats_df = pd.read_csv(stats_file)
    validation_percentage = stats_df["validation_percentage"].iloc[0]
    print(f"Loaded validation percentage for `{args.formula}`: {validation_percentage}%")
except FileNotFoundError:
    print(f"Validation stats file not found: {stats_file}")
    cursor.close()
    conn.close()
    exit(1)
except KeyError as e:
    print(f"Invalid file format, missing key: {e}")
    cursor.close()
    conn.close()
    exit(1)

# Update the `lipinski_valid_percent` column for the formula
update_percentage_query = """
UPDATE chemical_formulae 
SET lipinski_valid_percent = %s 
WHERE id = %s;
"""
try:
    cursor.execute(update_percentage_query, (validation_percentage, formula_id))
    print(f"Updated `lipinski_valid_percent` for `{args.formula}` to {validation_percentage}%.")
except pymysql.MySQLError as e:
    print(f"Error updating lipinski_valid_percent: {e}")

# Convert 'Fragments' column to valid JSON
def to_valid_json(fragment):
    try:
        return json.dumps(eval(fragment))
    except (SyntaxError, TypeError, NameError):
        return json.dumps([])

df["Fragments"] = df["Fragments"].apply(to_valid_json)

# Insert data into `lipinski_results` linked to the formula_id
for _, row in df.iterrows():
    insert_query = """
    INSERT INTO lipinski_results (formula_id, smiles, hba, hbd, logp, lipinski_valid, error, fragments)
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
    """
    try:
        cursor.execute(
            insert_query,
            (
                formula_id,
                row['SMILES'],
                row['HBA'],
                row['HBD'],
                row['LogP'],
                row['lipinski_valid'],
                row['Error'],
                row['Fragments'],
            )
        )
        print(f"Inserted row with SMILES: {row['SMILES']} under formula ID: {formula_id}")
    except pymysql.MySQLError as e:
        print(f"Error inserting row with SMILES {row['SMILES']}: {e}")

# Commit changes and close connection
try:
    conn.commit()
    print("All changes committed successfully.")
except pymysql.MySQLError as e:
    print(f"Error committing changes: {e}")

cursor.close()
conn.close()
print("Connection closed. Data insertion process complete.")
