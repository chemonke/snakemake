rule run_sql_container:
    input:
        config["files"]["sql_sentinel"],         # Ensure MariaDB image is pulled
        config["files"]["network_sentinel"]      # Ensure Docker network exists
    output:
        config["files"]["sql_running_sentinel"]
    params:
        root_password=config["sql"]["root_password"],
        database=config["sql"]["database"]
    shell:
        """
        docker stop mariadb || true
        docker rm mariadb || true
         # Allow root access from any host by setting MYSQL_ROOT_HOST=% 
        docker run --rm --name mariadb --network=my_network \
            -e MYSQL_ROOT_PASSWORD={params.root_password} \
            -e MYSQL_DATABASE={params.database} \
            -e MYSQL_ROOT_HOST=% \
            -d mariadb:10.5
        echo "Waiting for MariaDB to be ready..."
        for i in {{1..20}}; do
            docker exec mariadb mysqladmin ping -uroot -p{params.root_password} > /dev/null 2>&1 && break
            sleep 3
        done || (echo "MariaDB failed to start" && exit 1)
            # Additional check to ensure a query works
        for i in {{1..5}}; do
            if docker exec mariadb mysql -uroot -p{params.root_password} -e "SELECT 1;" > /dev/null 2>&1; then
                break
            fi
            sleep 3
        done || (echo "MariaDB did not respond to queries" && exit 1)
            touch {output}
        """


rule insert_data_to_db:
    input:
        f"{output}/{{formula}}_lipinski_results_with_fragments.csv",  # Input CSV file
        config["files"]["pydev_sentinel"],                             # Ensure pydev container is built
        config["files"]["sql_running_sentinel"]                       # Ensure MariaDB container is running
    output:
        f"{sentinel}/{{formula}}_data_inserted.flag"
    params:
        user="root",
        password=config["sql"]["root_password"],
        host="mariadb",  # Container name
        database=config["sql"]["database"]
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) --network=my_network -v $(pwd):/workspace -w /workspace pydev \
        /opt/conda/envs/pydev/bin/python scripts/sql.py \
        --input {input[0]} \
        --formula {wildcards.formula} \
        --user {params.user} \
        --password {params.password} \
        --host {params.host} \
        --database {params.database}
        touch {output}
        """
