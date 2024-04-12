import numpy as np
import os
import tarfile
from termcolor import colored

def bin2targz(directory, start_filter="kplr"):
    print("Searching in directory :" +  colored("{}".format(directory), "red"))
    # Step 1: List all subdirectores within "directory"
    if start_filter == None:
        subdirs = [subdir for subdir in os.listdir(directory) if os.path.isdir(os.path.join(directory, subdir))]
    else:
        subdirs = [subdir for subdir in os.listdir(directory) if os.path.isdir(os.path.join(directory, subdir)) and subdir.startswith(start_filter)]
    Ndirs=len(subdirs)
    print("Number of subdirectory found within : " + colored("{}".format(Ndirs), "red"))
    print("Processing each directory...")
    index=1
    for subdir in subdirs:
        print(" ----> [{}/{}] {}".format(index, Ndirs, subdir))
        subdir_path = os.path.join(directory, subdir)
        output_dir = os.path.join(subdir_path, "outputs")
        # Step 2: Compress each file as tar.gz
        for f in os.listdir(output_dir):
            if f.endswith(".bin"):
                targzfile=f.replace(".bin", ".tar.gz")
                str1=colored("      - File ", 'cyan')
                str2=colored("{}".format(f, targzfile), 'yellow')
                print(str1 + " ====> " + str2 + "...",  end='')
                file_path = os.path.join(output_dir, f)
                tar_file_path = os.path.join(output_dir, targzfile)
                with tarfile.open(tar_file_path, "w:gz") as tar:
                    tar.add(file_path, arcname=os.path.basename(file_path))
                # Step 3: Delete the original bin file
                os.remove(file_path)
                print(colored(" Compression done ", 'green') + "..." + colored(" Original bin file removed", 'red'))
        index=index+1


def bin2targz_conditional(directory, start_filter="kplr", phase="A"):
    print("Searching in directory :" +  colored("{}".format(directory), "red"), flush=True)
    # Step 1: List all subdirectores within "directory"
    if start_filter == None:
        subdirs = [subdir for subdir in os.listdir(directory) if os.path.isdir(os.path.join(directory, subdir))]
    else:
        subdirs = [subdir for subdir in os.listdir(directory) if os.path.isdir(os.path.join(directory, subdir)) and subdir.startswith(start_filter)]
    Ndirs=len(subdirs)
    print("Number of subdirectory found within : " + colored("{}".format(Ndirs), "red"), flush=True)
    print("Processing each directory...", flush=True)
    index=1
    for subdir in subdirs:
        process=False # Flag used to know whether we proceed to compression or not
        print(" ----> [{}/{}] {}".format(index, Ndirs, subdir), flush=True)
        subdir_path = os.path.join(directory, subdir)
        output_dir = os.path.join(subdir_path, "outputs")
        # Step 2: Look at the content of the params.hdr to get Nsamples and Nsamples_done
        params_files = [f for f in os.listdir(output_dir) if f.endswith(f"_{phase}_params.hdr")]
        if len(params_files) == 1:
            params_file = os.path.join(output_dir, params_files[0])
        else:
            print(colored("     Error: multiple params files found! This situation should not happen in real case", "red"), flush=True)
            print(colored("List of found files: ","red"), flush=True)
            for p in params_files:
                print("    - ", p, flush=True)
            print(colored("     Debug required. The program will stop now", "red"), flush=True)
            exit()
        if os.path.isfile(params_file):
            with open(params_file, "r") as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith("! Nsamples_done="):
                        nsamples_done = int(line.split("=")[1].strip())
                    if line.startswith("! Nsamples="):
                        Nsamples=int(line.split("=")[1].strip())
            if nsamples_done == Nsamples - 1:
                print(colored("     hdr file found with Nsamples = nsamples_done = {}".format(Nsamples), "green"), flush=True)
                process=True # All conditions met to compress the bin files into tar.gz
            else:
                print(colored("     hdr file found but does not met the criteria: nsamples_done ={} != {}".format(nsamples_done+1, Nsamples), "yellow"), flush=True)
        else:
            print(colored("     No file of type: ", f"*_{phase}_params.hdr in this directory","yellow"), flush=True)    

        # Step 3: Compress each file as tar.gz if the conditions are met
        if process == True:
            for f in os.listdir(output_dir):
                if f.endswith(".bin"):
                    targzfile=f.replace(".bin", ".tar.gz")
                    str1=colored("      - File ", 'cyan')
                    str2=colored("{}".format(f, targzfile), 'yellow')
                    print(str1 + " ====> " + str2 + "...",  end='', flush=True)
                    file_path = os.path.join(output_dir, f)
                    tar_file_path = os.path.join(output_dir, targzfile)
                    with tarfile.open(tar_file_path, "w:gz") as tar:
                        tar.add(file_path, arcname=os.path.basename(file_path))
                    # Step 3: Delete the original bin file
                    os.remove(file_path)
                    print(colored(" Compression done ", 'green') + "..." + colored(" Original bin file removed", 'red'), flush=True)
        index=index+1

def test_bin_realdata(dir_root="../data/outputs/1001/", enforce_nsamples_is_nsamples_done=True):
    # Usage example of the two implementations
    if enforce_nsamples_is_nsamples_done == False:
        print("Testing the bin2targz function...")
        bin2targz(dir_root, start_filter="kplr") # This compress any file of .bin extension that match the filter
    else:
        print("Testing the bin2targz_conditional     function...")
        bin2targz_conditional(dir_root, start_filter="kplr", phase="A") # This compress ONLY if Nsamples = Nsamples_done in the *_params.hdr