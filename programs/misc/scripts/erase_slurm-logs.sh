#!/bin/bash

# Deletes all files that are of syntax "slurm-*.out" within a whole tree of directories

# Check if a root directory argument is provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <root_directory>"
  exit 1
fi

# Get the root directory from the command line argument
root_directory="$1"

# Use the find command to search for files matching the specified criteria
find "$root_directory" -type f -name "slurm-*.out" -delete

