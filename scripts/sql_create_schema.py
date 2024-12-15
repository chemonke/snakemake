import pymysql
import argparse
import os
print(f"Executing schema creation from: {os.path.abspath(__file__)}")


# Parse arguments for database credentials
parser = argparse.ArgumentParser(description="Create database schema")
parser.add_argument("--host", required=True, help="Database host")
parser.add_argument("--user", required=True, help="Database username")
parser.add_argument("--password", required=True, help="Database password")
parser.add_argument("--database", required=True, help="Database name")
args = parser.parse_args()

# Database connection settings
db_config = {
    "host": args.host,
    "user": args.user,
    "password": args.password,
    "database": args.database,
}

# Table creation queries
TABLE_CREATION_QUERIES = [
    """
    CREATE TABLE IF NOT EXISTS chemical_formulae (
        id INT AUTO_INCREMENT PRIMARY KEY,
        formula TEXT NOT NULL,
        mw FLOAT DEFAULT NULL,
        description TEXT DEFAULT NULL,
        valid_percent FLOAT DEFAULT NULL,
        UNIQUE KEY (formula(255))
    ) ENGINE=InnoDB;
    """,
    """
    CREATE TABLE IF NOT EXISTS computational_results (
        id INT AUTO_INCREMENT PRIMARY KEY,
        formula_id INT NOT NULL,
        smiles TEXT NOT NULL UNIQUE,
        hba INT,
        hbd INT,
        logp FLOAT,
        valid BOOLEAN,
        lipinski_valid BOOLEAN,
        error BOOLEAN,
        fragments TEXT,
        chemval BOOLEAN,  -- Add this column
        chem_violation TEXT,  -- Add this column
        FOREIGN KEY (formula_id) REFERENCES chemical_formulae(id) ON DELETE CASCADE
    ) ENGINE=InnoDB;   
    """,
    """
    CREATE TABLE IF NOT EXISTS web_data (
        id INT AUTO_INCREMENT PRIMARY KEY,
        computational_results_id INT NOT NULL,
        cas_number TEXT NOT NULL UNIQUE,
        vendor_available BOOLEAN,
        patent_available BOOLEAN,
        retrieved_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (computational_results_id) REFERENCES computational_results(id) ON DELETE CASCADE
    ) ENGINE=InnoDB;
    """,
    """
    CREATE TABLE IF NOT EXISTS vendors (
        id INT AUTO_INCREMENT PRIMARY KEY,
        web_data_id INT NOT NULL,
        vendor_name TEXT,
        price FLOAT,
        purity FLOAT,
        availability INT,
        FOREIGN KEY (web_data_id) REFERENCES web_data(id) ON DELETE CASCADE
    ) ENGINE=InnoDB;
    """
]

# Connect and create tables
try:
    conn = pymysql.connect(**db_config)  # Use db_config here
    cursor = conn.cursor()
    for query in TABLE_CREATION_QUERIES:
        print(f"Executing query:\n{query}")  # Log each query
        cursor.execute(query)
    print("Database schema created successfully.")
except pymysql.MySQLError as e:
    print(f"Error creating database schema: {e}")
finally:
    cursor.close()
    conn.close()
