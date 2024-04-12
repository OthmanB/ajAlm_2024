#!/bin/bash

# Get the directory to search for subdirectories from the command line argument
directory="$1"

# Loop through all subdirectories within the directory
for subdir in "$directory"/*; do
    if [ -d "$subdir" ]; then
        # Erase files of syntax "*_B_*.*", excluding _restore_ files and files under 'restore' directory
        echo "Removing B files in subdirectory: ${subdir##*/}"
        find "$subdir" -type f -name "*_B_*.*" ! -name "*_restore_*" ! -path "*/restore/*" -delete 

        # Find and compress the entire subdirectory if it contains "L" files
        if [ -n "$(find "$subdir" -type f -name "*_L_*.*" -print -quit)" ]; then
            echo "Compressing subdirectory: ${subdir##*/}"
            tar -czf "$subdir.tar.gz" -C "$directory" "${subdir##*/}"
            # Remove the "L" files within the subdirectory's outputs directory
            echo "Removing 'L' files in subdirectory: ${subdir##*/}"
            find "$subdir/outputs/" -type f -name "*_L_*.bin" -delete
        fi
    fi
done

