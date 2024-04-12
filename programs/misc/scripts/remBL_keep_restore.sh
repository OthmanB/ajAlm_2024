#!/bin/bash

# Get the directory to search for subdirectories from command line argument
directory="$1"

# Loop through all subdirectories within the directory
for subdir in "$directory"/*; do
    if [ -d "$subdir" ]; then
        # Erase files of syntax "*_B_*.*", excluding _restore_ files and files under 'restore' directory
        echo "Removing B files in subdirectory: ${subdir##*/}"
        find "$subdir" -type f -name "*_B_*.*" ! -name "*_restore_*" ! -path "*/restore/*" -delete 

        # Compress the files of syntax "*_L_*.*" in a tar.gz file
        #echo "Compressing L files in subdirectory: ${subdir##*/}"
        #find "$subdir/" -name "*_L_*.*" -print0 | tar -czvf "$subdir.tar.gz" --null -T -

        # Remove the files of syntax "*_L_*.*" after compressing, excluding _restore_ files and files under 'restore' directory
        echo "Removing L files in subdirectory: ${subdir##*/}"
        find "$subdir" -type f -name "*_L_*.*" ! -name "*_restore_*" ! -path "*/restore/*" -delete
    fi
done

