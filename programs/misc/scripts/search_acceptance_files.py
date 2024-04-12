import os
from termcolor import colored
import fnmatch

def show_colored_line(line, range_min, range_max, i0=0, iend=None, separators=(' ', '\t'), prefix=""):
    elements = [element for element in line.split(separators[0]) if element]
    if len(elements) == 1:
        elements = [element for element in line.split(separators[1]) if element]
    parsed_elements = [float(element) for element in elements]
    if iend is None:
        iend = len(parsed_elements) - 1

    for i in range(len(parsed_elements)):
        if i < i0 or i > iend:
            print(prefix + colored(elements[i], 'white'), end=' ')
        elif parsed_elements[i] < range_min or parsed_elements[i] > range_max:
            print(prefix + colored(elements[i], 'red'), end=' ')
        else:
            print(prefix + colored(elements[i], 'green'), end=' ')
    print()  # Print a new line

    for i in range(len(parsed_elements)):
        if i < i0 or i > iend:
            return False
        elif parsed_elements[i] < range_min or parsed_elements[i] > range_max:
            return False
    return True


def read_pattern_files(directory, file_pattern="*"):
    '''
        Scan a given directory and its subdirectories for files that match the specified phase pattern.
        Returns a list of tuples that contain (filename, file_path, last_line).
        The last line of an acceptance file contains the last acceptance rate for all parallel chains
    '''
    found_files = []
    info_all = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if fnmatch.fnmatch(file, file_pattern):
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as f:
                    info = {}
                    lines = f.readlines()
                    if lines:
                        last_line = lines[-1].strip()
                        info["last_line"] = last_line
                        found_files.append((file, file_path))
                        info_all.append(info)
    return found_files, info_all

def display_file_paths(files_list, info=None, highlight_params=[1, None, 0.15, 0.35]):
    # Create a dictionary to store the branches
    branches = {}
    max_shift_spaces = 0  # Variable to track the maximum shift spaces needed
    for i, file_tuple in enumerate(files_list):
        filename, file_path = file_tuple
        # Normalize file path for consistency
        file_path = os.path.normpath(file_path)
        # Split file_path into directory components
        directories = file_path.split(os.sep)
        # Get the key for the branch
        branch_key = os.sep.join(directories[:-1])
        # Get the branch from the dictionary, or create a new one if it doesn't exist
        if branch_key in branches:
            branch = branches[branch_key]
        else:
            branch = []
            branches[branch_key] = branch
        # Add the file to the branch
        if info is not None:
            branch.append((directories[-1], info[i]))
            shift_spaces = len(directories[-1]) - len(filename) - 4
            max_shift_spaces = max(max_shift_spaces, shift_spaces)
        else:
            branch.append((directories[-1]))
    # Traverse the branches and display the tree
    for branch_key, branch in branches.items():
        # Split the branch_key into directories
        directories = branch_key.split(os.sep)
        num_directories = len(directories)
        # Print the graphical representation of the branch
        for i, directory in enumerate(directories):
            if i == num_directories - 1:
                # Last directory
                print(f"{directory}/")
            elif i == 0:
                # First directory
                print(f"+--{directory}/")
            else:
                # Intermediate directories
                print(f"|  " * i + f"+--{directory}/")
        # Print the files in the branch
        if info is not None:
            shift_spaces=num_directories
            for file_name, info_name in branch:
                print(f"|  " * num_directories + f"+--{file_name}")
                shift_spaces = len(file_name) - len("slurm") + max_shift_spaces
                line = info_name["last_line"]
                line_highlighted=show_colored_line(line, highlight_params[2], highlight_params[3], 
                                                   i0=highlight_params[0], 
                                                   iend=highlight_params[1],
                                                   prefix="     ")
                                                   #prefix=" " * shift_spaces)
        else:
            for file_name in branch:
                print(f"|  " * num_directories + f"+--{file_name}")



def show_acceptance_files(directory, phase="A"):
    print("Searching for files with pattern : ", "*_" + phase + "_acceptance.txt")
    found_files, infos=read_pattern_files(directory, file_pattern="*_" + phase + "_acceptance.txt")
    if found_files == []:
        print("We could not find any file with that syntax!")
    else:
        display_file_paths(found_files, info=infos)


#root_dir="/scratch/ob19/ajAlm_project/data/TRANSFORMED_RUN/Outputs/Legacy/"
root_dir="/scratch/ob19/ajAlm_project/data/TRANSFORMED_RUN/Outputs/Kamiaka2018/"
show_acceptance_files(root_dir, phase="A")
