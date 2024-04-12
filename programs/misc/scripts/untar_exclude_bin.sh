#!/bin/bash

# This script allows you to untar the content an ensemble of tar files
# but excluding the files with a bin extension

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <source_directory> <destination_directory>"
    exit 1
fi

# Assign command-line arguments to variables
source_directory="$1"
destination_directory="$2"

# Ensure the source directory exists
if [ ! -d "$source_directory" ]; then
    echo "Error: Source directory '$source_directory' not found."
    exit 1
fi

# Ensure the destination directory exists, or create it
if [ ! -d "$destination_directory" ]; then
    mkdir -p "$destination_directory"
fi

# Iterate over all tar.gz files in the source directory
for tar_file in "$source_directory"/*.tar.gz; do
    # Extract files, excluding those with .bin extension, into the destination directory
    echo "Extracting $tar_file to $destination_directory"
    tar --exclude='*.bin' -xzvf "$tar_file" -C "$destination_directory"
done

