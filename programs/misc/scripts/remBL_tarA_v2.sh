#!/bin/bash

# Example of use: 
#    Delete B, L, compress remaining (A) and delete directory after :  ./your_script.sh -d /path/to/directory
#    Delete B, L, compress remaining (A) BUT DO NOT DELETE directory: ./your_script.sh /path/to/directory

# Initialize variables
directory=""
delete_directory=false

# Parse command line options
while getopts ":d" opt; do
    case $opt in
        d)
            delete_directory=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Shift the command line options so that "$1" now refers to the directory argument
shift $((OPTIND-1))

# Get the directory to search for subdirectories from command line argument
directory="$1"

# Loop through all subdirectories within the directory
for subdir in "$directory"/*; do
    if [ -d "$subdir" ]; then
        # Erase files of syntax "*_L_*.*" and "*_B_*.*"
        echo "Removing files in subdirectory: ${subdir##*/}"
        find "$subdir" -type f -name "*_L_*.*" -o -name "*_B_*.*" -delete        
        # Compress the subdirectory in a tar.gz file
        echo "Compressing subdirectory: ${subdir##*/}"
        tar -czf "$subdir.tar.gz" -C "$directory" "${subdir##*/}"
        
        if [ "$delete_directory" = true ]; then
            # Remove the subdirectory if the option is enabled
            echo "Removing subdirectory: ${subdir##*/}"
            rm -r "$subdir"
        fi
    fi
done

echo "All Done"
