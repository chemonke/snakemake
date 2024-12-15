# debugging option; clean outputs and sentinel files for re-running the workflow
# note that debugging must be enabled in the config (safety feature to prevent accidental deletion)

rule clean:
    shell:
        """
        if [ "{config[debug]}" = "True" ]; then
            echo "Cleaning up outputs and sentinel files for debugging..."
            rm -rf {config[folders][output]}/* {config[folders][sentinel]}/*
        else
            echo "Debug mode is off. Cleanup skipped."
        fi
        """
