rule build_pydev_image:
    output:
        config["files"]["pydev_sentinel"]
    shell:
        """
        if ! docker images | grep -q 'pydev'; then
            docker build -t pydev {containers}/python-container
        fi
        touch {output}
        """

rule build_surge_image:
    output:
        config["files"]["surge_sentinel"]
    shell:
        """
        if ! docker images | grep -q 'surge'; then
            docker build -t surge {containers}/surge
        fi
        touch {output}
        """

rule build_sql_image:
    output:
        config["files"]["sql_sentinel"]
    shell:
        """
        if ! docker images | grep -q 'mariadb'; then
            docker pull mariadb:10.5  # Replace 10.5 with the desired MariaDB version
        fi
        touch {output}
        """

rule create_docker_network:
    output:
        config["files"]["network_sentinel"]
    shell:
        """
        if ! docker network ls | grep -q 'my_network'; then
            docker network create my_network
        fi
        touch {output}
        """
