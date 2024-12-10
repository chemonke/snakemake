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
