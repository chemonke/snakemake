from flask import Flask, render_template
import pymysql

app = Flask(__name__)

# Database connection settings
DB_CONFIG = {
    "host": "mariadb",  # Use the container's hostname
    "user": "root",
    "password": "swoosh",
    "database": "leads",
}


def check_database_status():
    """
    Check if the database is running and accessible.
    Returns True if connected, otherwise False.
    """
    try:
        conn = pymysql.connect(**DB_CONFIG)
        conn.close()
        return True
    except pymysql.MySQLError:
        return False


@app.route("/")
def index():
    # Check the database status
    db_status = check_database_status()
    return render_template("index.html", db_status=db_status)


if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=5000)

