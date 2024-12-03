#!/usr/bin/env bash

# Usage: ./run_surge.sh input_string output_file
INPUT_STRING=$1
OUTPUT_FILE=$2

docker run --rm -v "$(pwd)":/data --user $(id -u):$(id -g) surge:latest surge -S -o/data/"$OUTPUT_FILE" "$INPUT_STRING"
