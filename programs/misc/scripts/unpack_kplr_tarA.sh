#!/bin/bash

# Example of use:
#    Decompress tar and DELETE the tar.gz after :  ./your_script.sh -d /path/to/directory
#    Decompress tar BUT DO NOT DELETE the tar.gz file : ./your_script.sh /path/to/directory

# Initialize variables
directory=""
delete_tarfile=false

# Parse command line options
while getopts ":d" opt; do
    case $opt in
        d)
            delete_tarfile=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Shift the command line options so that "$1" now refers to the directory argument
shift $((OPTIND-1))

# Get the directory to search for tar.gz files from the command line argument
directory="$1"

# Function to unpack tar.gz files in a directory and its subdirectories
unpack_tar_files() {
    local dir="$1"
    for tarfile in "$dir"/*.tar.gz; do
        if [ -f "$tarfile" ] && [[ $(basename "$tarfile") == kplr* ]]; then
            echo "Unpacking $tarfile in $dir"
            tar -xzf "$tarfile" -C "$dir"
            if [ "$delete_tarfile" = true ]; then
                echo "Deleting $tarfile"
                rm "$tarfile"
            fi
        fi
    done

    for subdir in "$dir"/*; do
        if [ -d "$subdir" ]; then
            unpack_tar_files "$subdir"
        fi
    done
}

# Start the unpacking process in the specified directory
unpack_tar_files "$directory"

