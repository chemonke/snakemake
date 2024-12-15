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
        validation_stats=f"{output}/{{formula}}_validation_error_stats.csv",
        validation_plot=f"{output}/{{formula}}_validation_error_visualization.png",
        fragment_analysis=f"{output}/{{formula}}_fragment_analysis.csv",
        fragment_plot=f"{output}/{{formula}}_fragment_frequency_plot.png",
        molecules=f"{output}/{{formula}}_molecules.png",
        summary_stats=f"{output}/{{formula}}_summary_statistics.csv",
        pairplot=f"{output}/{{formula}}_pairplot.png",
        results=f"{output}/{{formula}}_results.csv",
        sql_sentinel=config["files"]["sql_running_sentinel"],
        pydev_sentinel=config["files"]["pydev_sentinel"]
    output:
        f"{sentinel}/{{formula}}_data_inserted.flag"
    params:
        user="root",  # Add your database user here
        password=config["sql"]["root_password"],  # Reference from config.yaml
        host="mariadb",  # Docker container hostname
        database=config["sql"]["database"]  # Reference from config.yaml
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) --network=my_network \
            -v $(pwd):/workspace -w /workspace pydev \
            /opt/conda/envs/pydev/bin/python scripts/sql.py \
            --input {input.results} \
            --stats {input.validation_stats} \
            --formula {wildcards.formula} \
            --user {params.user} \
            --password {params.password} \
            --host {params.host} \
            --database {params.database}
        touch {output}
        """
