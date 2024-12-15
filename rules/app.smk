rule app:
    params:
        port="5000",  # Port to bind Flask app
    shell:
        """
        docker run --rm --user $(id -u):$(id -g) \
            -v $(pwd)/app:/workspace/app -v $(pwd)/app/templates:/workspace/app/templates \
            -w /workspace/app -p {params.port}:{params.port} pydev \
            /opt/conda/envs/pydev/bin/python app.py
        """
