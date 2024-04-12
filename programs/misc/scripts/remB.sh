#!/bin/bash

# Get the directory to search for subdirectories from command line argument
directory="$1"

# Loop through all subdirectories within the directory
for subdir in "$directory"/*; do
    if [ -d "$subdir" ]; then
        # Erase files of syntax "*_B_*.*"
        echo "Removing B files in subdirectory: ${subdir##*/}"
        find "$subdir" -type f -name "*_B_*.*" -delete
    fi
done
