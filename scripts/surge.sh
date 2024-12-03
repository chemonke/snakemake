#!/usr/bin/env bash

# Usage: ./run_surge.sh input_string output_file
INPUT_STRING=$1
OUTPUT_FILE=$2

# Create a named pipe
PIPE_NAME=$(mktemp -u)
mkfifo "$PIPE_NAME"

SAMPLE_SIZE=10000

# Run `surge` inside a `timeout` command to terminate after 5 minutes
timeout 5m docker run --rm -v "$(pwd)":/data --user $(id -u):$(id -g) surge:latest surge -S "$INPUT_STRING" > "$PIPE_NAME" &

# Capture the PID of the background process
DOCKER_PID=$!

# Use `awk` to perform random sampling and limit output to n lines
awk 'BEGIN {srand()} {if (rand() < 0.01) print $0}' < "$PIPE_NAME" | head -n $SAMPLE_SIZE > "$OUTPUT_FILE"

# Wait for the `docker run` process to complete, but handle the timeout case
wait $DOCKER_PID

# Check if `timeout` terminated the process
if [[ $? -eq 124 ]]; then
  echo "Process terminated after exceeding 5 minutes."
fi

# Clean up the named pipe
rm "$PIPE_NAME"
