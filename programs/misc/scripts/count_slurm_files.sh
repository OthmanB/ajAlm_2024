#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <root_directory>"
    exit 1
fi

root_directory="$1"

# Find all directories containing a file named cpptamcmc under the specified root directory
directories=$(find "$root_directory" -type f -name cpptamcmc -exec dirname {} \;)

# Iterate over the found directories
for dir in $directories; do
    # Count the number of files starting with "slurm" in the current directory
    count=$(find "$dir" -type f -name "slurm*" | wc -l)
    
    # Print the directory path and the file count
    echo "Directory: $dir"
    echo "Number of files starting with 'slurm': $count"
done

