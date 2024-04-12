#!/bin/bash

# Directory where you want to search for and delete the files
search_directory="/scratch/ob19/ajAlm_project/data/TRANSFORMED_RUN/Outputs/Kamiaka2018/"

# Use the find command to locate and delete the files
find "$search_directory" -type f -name "kplr*.tar.gz" -exec rm -f {} \;

echo "Files with 'kplr' prefix and ending in '.tar.gz' deleted."

