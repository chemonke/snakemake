
import pymysql
import pandas as pd
import json  # For JSON validation
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Insert data into MariaDB")
parser.add_argument("--input", required=True, help="Path to the input CSV file")
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

# Create table if not exists
# Create table if not exists
create_table_query = """
CREATE TABLE IF NOT EXISTS lipinski_results (
    id INT AUTO_INCREMENT PRIMARY KEY,
    smiles TEXT NOT NULL,
    mw FLOAT,
    hba INT,
    hbd INT,
    logp FLOAT,
    valid BOOLEAN,
    error BOOLEAN,
    fragments JSON,
    UNIQUE KEY (smiles(255))  -- Enforces uniqueness for smiles column
);
"""

try:
    cursor.execute(create_table_query)
    print("Table `lipinski_results` is ready.")
except pymysql.MySQLError as e:
    print(f"Error creating table: {e}")
    cursor.close()
    conn.close()
    exit(1)

# Load the CSV into a DataFrame
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

# Validate required columns
required_columns = {'SMILES', 'MW', 'HBA', 'HBD', 'LogP', 'Valid', 'Error', 'Fragments'}
if not required_columns.issubset(df.columns):
    missing_columns = required_columns - set(df.columns)
    print(f"Missing required columns: {missing_columns}")
    cursor.close()
    conn.close()
    exit(1)

# Convert 'Fragments' column to valid JSON
def to_valid_json(fragment):
    try:
        # Convert string to Python list and then to JSON string
        return json.dumps(eval(fragment))
    except (SyntaxError, TypeError, NameError):
        # Handle invalid formats gracefully
        return json.dumps([])

df["Fragments"] = df["Fragments"].apply(to_valid_json)

# Insert data into the table
for _, row in df.iterrows():
    insert_query = """
    INSERT INTO lipinski_results (smiles, mw, hba, hbd, logp, valid, error, fragments)
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
    """
    try:
        cursor.execute(
            insert_query,
            (
                row['SMILES'],
                row['MW'],
                row['HBA'],
                row['HBD'],
                row['LogP'],
                row['Valid'],
                row['Error'],
                row['Fragments'],
            )
        )
        print(f"Inserted row with SMILES: {row['SMILES']}")
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