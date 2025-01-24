#!/usr/bin/env bash

# Define the input path
TARGET_IN=$1

# Define the output file
TARGET_OUT=$2

# Clear the output file if it already exists
> "$TARGET_OUT"

# Loop through all .log files in the current directory
for file in "$TARGET_IN"/*.dynamics; do
    if [ -f "$file" ]; then
        echo "$file"
        tail +2 "$file" >> "$TARGET_OUT"
    fi
done
