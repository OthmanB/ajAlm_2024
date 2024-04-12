import os
import re

import numpy as np

def find_max_prior(filename):
    '''
    This function is designed to read *.priors files and get the max of 
    the joint-PDF in the (x,y) plane 
    '''
    # Read the table from the file
    with open(filename, 'r') as f:
        table = f.read()

    # Ignore lines starting with "#" or empty lines
    lines = table.split('\n')
    lines = [line.strip() for line in lines if line and not (line.startswith("#") or line.startswith("!") or line.startswith("*"))]
    
    # Extract x-axis values
    x_axis = np.array(lines[0].split()[1:])
    #print("x_axis :", x_axis)

    # Extract y-axis values and joint values
    y_axis = []
    joint_values = []
    for line in lines[1:]:
        values = line.split()
        #print("values : ", values)
        y_axis.append(values[0])
        joint_values.append(list(map(float, values[1:])))
    joint_values = np.array(joint_values)
    #print("y_axis : ", y_axis)

    # Find the maximum joint value and its corresponding indices
    max_joint_value = np.max(joint_values)
    max_joint_indices = np.where(joint_values == max_joint_value)
    max_joint_indices = (max_joint_indices[0][0], max_joint_indices[1][0])
    print(max_joint_indices)
    # Get the corresponding (xmax, ymax) values
    xmax = x_axis[max_joint_indices[1]]
    ymax = y_axis[max_joint_indices[0]]
    
    return xmax, ymax, max_joint_value



def read_a2a3a4_prior_table(filename):
    '''
        Priors for a2,a4 are calculated by the other program in
        programs/acoefs2activity/eval_aj_range.py 
        and a3 coeficient are set to some threshold depending on a1
        This is saved in a table that can be read here and used to reformat
        inputs model files
    '''
    header = []  # List to store the commented lines
    lines = []   # List to store the non-commented lines
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line == '':
                header.append(line)
            else:
                row = line.split()
                lines.append(row)
    return header, lines

def get_line_with_key(table, key, i=0):
    '''
        A small function to retrieve the input that match a key at index i
        Used to search within tables read by read_a2a3a4_prior_table()
    '''
    for line in table:
        #print(" line : ", line)
        #print("        line[{}] = {}".format(i, line[i]))
        if line[i] == key:
            return line
    return None

def find_modelfile_byID(directory, identifier, prefix="kplr", suffix="", extension=".model"):
    '''
        Search for a model file that within a directory, provided its numeric identifier (eg. KIC)
        and assuming a prefix and suffix
    '''
    for filename in os.listdir(directory):
        if filename.startswith(prefix) and filename.endswith(suffix+extension) and identifier in filename:
            return filename
    raise FileNotFoundError("No file found with the given syntax.")

def format_string_length(s, len_p, req_len_p, padding_char=' '):
    '''
        Function that pads a string to a desired length
    '''
    if req_len_p < len_p:
        raise ValueError("The desired length is smaller than the current length.")
    formatted_string = s.rjust(req_len_p, padding_char)
    return formatted_string

def search_index(lines, key, stop_not_found=False, show_warning=True):
    '''
        Small code that look for a location 'key' (eg. "freq_smoothness") inside 
        an ensemble of lines in an array
    '''
    index=None
    for i, line in enumerate(lines[:]):
        words= line.split() # Look for only words separated by spaces or tab... This avoid True when eg. model_MS_Global_a1etaa3_HarveyLike
        if key in words:
            index = i
            break  
    if index == None:
        if show_warning == True:
            print("Warning: Key '{}' not found".format(key))
    if stop_not_found == True:
        raise ValueError("Error : Key '{}' not found in provided lines. Debug required.".format(key))
    
    return index

def polish_file(file_path, is_ajAlm_model=False):
    '''
        This function allows to update the outdated .model files that were used until 2020
        into the new format, compatible with TAMCMC version >1.80.0 
        This function only ensure that we have the standard and most up to date set of global
        priors in the file. 
        This function DOES NOT save the file. It returns a structured template that has to be written by
        another function. 
        Because it creates a compatible template, this function does not change the inputs and priors. 
        So you will likely need to combine it with another code to change the 
        input values, as it sets a default set that is of Fix type when the variable was not found.
        Thus, Please use polish_file() in conjonction with update_file()
    '''
    # Presets
    if is_ajAlm_model == False:
        modelname = "model_MS_Global_aj_HarveyLike"
    else:
        modelname = "model_MS_Global_ajAlm_HarveyLike"
        decompose_Alm="1"
        filter_type="gate"

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lines0=lines
    # Find the index of the label
    label_index = None
    for i, line in enumerate(lines[:]):
        if line.startswith("# Controls and priors for common parameters"):
            label_index = i
            break

    # Remove 'Asphericity_eta' lines
    #lines = [line for line in lines if "Asphericity_eta" not in line.split()]
    lines = [line for line in lines if "Asphericity_eta" not in line]

    if label_index is None:
        # Label not found, raise an error
        print(" Error while processing file :", file_path)
        raise(ValueError, "Could not found the line '# Controls and priors for common parameters'")
    
    # Add new model name, if required
    modelname_exists=False
    for i, line in enumerate(lines[:]):
        if "model_fullname" in line:
            modelname_exists=True
            #if is_ajAlm_model == False:
            elements = line.split()
            elements[-1]=modelname
            lines[i] = "            ".join(elements) + "\n"

    if modelname_exists == False:
        lines.insert(label_index + 1, "           model_fullname             {}\n".format(modelname))
    # Rename 'Splitting_a1' to 'a1_0' and a3, to be sure that we have similarly formated keywords
    for i, line in enumerate(lines[:]):
        words= line.split() # Look for only words separated by spaces or tab... This avoid True when eg. model_MS_Global_a1etaa3_HarveyLike
        if "Splitting_a1" in words or "a1" in words:
            #print("lines[{}] = {}".format(i, lines[i]) )
            lines[i] = line.replace("Splitting_a1", "        a1_0")
            lines.insert(i+1,"                     a1_1             Fix              0.000000\n")
            #print("lines[{}] = {}".format(i, lines[i]) )
            break
    for i, line in enumerate(lines[:]): # We iterate over a copy of lines here to avoid infinite loops
        words= line.split() # Look for only words separated by spaces or tab... This avoid True when eg. model_MS_Global_a1etaa3_HarveyLike
        if "Splitting_a3" in words or "a3" in words:
            lines[i] = line.replace("Splitting_a3", "        a3_0")
            lines.insert(i+1,"                     a3_1             Fix              0.000000\n")
            break
    # Check what is there and what is not
    a1_exists=0
    a2_exists=0
    a3_exists=0
    a4_exists=0
    a5_exists=0
    a6_exists=0
    epsilon_exists=0
    theta0_exists=False
    delta_exists=False
    decompose_Alm_exists=False
    filter_type_exists=False
    for i, line in enumerate(lines):
        words= line.split()
        if "a1_0" in words:
            a1_exists=a1_exists + 1
        if "a1_1" in words:
            a1_exists=a1_exists + 1
        if "a2_0" in words:
            a2_exists=a2_exists + 1
        if "a2_1" in words:
            a2_exists=a2_exists + 1
        if "a3_0" in words:
            a3_exists=a3_exists + 1
        if "a3_1" in words:
            a3_exists=a3_exists + 1
        if "a4_0" in words:
            a4_exists=a4_exists + 1
        if "a4_1" in words:
            a4_exists=a4_exists + 1
        if "a5_0" in words:
            a5_exists=a5_exists + 1
        if "a5_1" in words:
            a5_exists=a5_exists + 1
        if "a6_0" in words:
            a6_exists=a6_exists + 1
        if "a6_1" in words:
            a6_exists=a6_exists + 1
        # Specific to ajAlm models
        if "epsilon_0" in words:
            epsilon_exists=epsilon_exists + 1
        if "epsilon_1" in words:
            epsilon_exists=epsilon_exists + 1   
        if "theta0" in words:
            theta0_exists=True 
        if "delta" in words:
            delta_exists=True   
        if "decompose_Alm" in words:
            decompose_Alm_exists=True   
        if "filter_type" in words:
            filter_type_exists=True 
    # Verify that we do have a1_0 and a3_0 now... In case the file had nothing about that, we need to do so
    if a1_exists != 2:
        raise ValueError("Error: Could not find the correct number of a1 parameters in the source file. Manual debug required.")
    if a3_exists != 2:
        raise ValueError("Error: Could not find the correct number of a3 parameters source file. Manual debug required.")

    # Insert new lines at the appropriate index for a2 and a4
    if modelname == "model_MS_Global_aj_HarveyLike":
        # Then we just need to add a2 after a1_1, if not there
        index=search_index(lines, "a1_1", stop_not_found=False)        
        new_lines=[]
        if a2_exists == 0:
            new_lines.append("                     a2_0             Fix              0.000000\n")
            new_lines.append("                     a2_1             Fix              0.000000\n")
        if a2_exists == 1 or a2_exists > 2:
            raise ValueError("Error: Expected to find 0 or 2 arguments for a2. But found 1 or >2")
        lines = lines[:index+1] + new_lines + lines[index+1:]   
        index=search_index(lines, "a3_1", stop_not_found=False)        
        new_lines=[]
        if a4_exists == 0:
            new_lines.append("                     a4_0             Fix              0.000000\n")
            new_lines.append("                     a4_1             Fix              0.000000\n")
        if a4_exists == 1 or a4_exists > 2:
            raise ValueError("Error: Expected to find 0 or 2 arguments for a2. But found 1 or >2")
        lines = lines[:index+1] + new_lines + lines[index+1:]   
    # Insert new lines at the appropriate index for a5 and a6
    if a5_exists == 0:
        if modelname == "model_MS_Global_aj_HarveyLike":
            index=search_index(lines, "a4_1", stop_not_found=False)  
        else:
            index=search_index(lines, "a3_1", stop_not_found=False)
        new_lines = ["                     a5_0             Fix              0.000000\n",
                     "                     a5_1             Fix              0.000000\n"]
        lines = lines[:index+1] + new_lines + lines[index+1:]
    if a5_exists == 1 or a5_exists > 2:
        raise ValueError("Error: Expected to find 0 or 2 arguments for a5. But found 1 or >2")
    
    if a6_exists == 0:
        index=search_index(lines, "a5_1", stop_not_found=False)        
        new_lines = ["                     a6_0             Fix              0.000000\n",
                     "                     a6_1             Fix              0.000000\n"]
        lines = lines[:index+1] + new_lines + lines[index+1:]
    if a6_exists == 1 or a6_exists > 2:
        raise ValueError("Error: Expected to find 0 or 2 arguments for a6. But found 1 or >2")
    # Modify 'Asymetry' line
    asym_exist=False
    for i, line in enumerate(lines[:]):
        words= line.split() # Look for only words separated by spaces or tab... This avoid True when eg. model_MS_Global_a1etaa3_HarveyLike
        if "Asymetry" in words:
            lines[i] = "                 Asymetry             Jeffreys_abs     10.000000          5.000000          200.0000\n"
            asym_exist=True
            break
    if asym_exist == False: # If no asymetry was there, add it before "Inclination"
        for i, line in enumerate(lines[:]):
            words= line.split() # Look for only words separated by spaces or tab... This avoid True when eg. model_MS_Global_a1etaa3_HarveyLike
            if "Inclination" in words:
                lines.insert(i-1, "                 Asymetry             Jeffreys_abs     10.000000          5.000000          200.0000\n")    
    
    # Case of a2 or not a2 (model_MS_Global_aj_HarveyLike or model_MS_Global_ajAlm_HarveyLike)
    # This is made at the end because we need to be sure that a6_1 entry is here to insert the ajAlm case         
    if modelname == "model_MS_Global_ajAlm_HarveyLike":
        # Then we need to add epsilon_nl, theta0, delta, decompose_Alm and filter_type
        # Locate Insertion point for filter_type and decompose_Alm and insert them
        new_lines=[]
        if decompose_Alm_exists == False:
            new_lines.append("            decompose_Alm             Fix              {}\n".format(decompose_Alm))
        if filter_type_exists == False:
            new_lines.append("              filter_type             {}\n".format(filter_type))
        index=search_index(lines, "freq_smoothness", stop_not_found=False)           
        lines = lines[:index+1] + new_lines + lines[index+1:]
        # Locate Insertion point for epsilon, theta and epsilon and insert them if required
        index=search_index(lines, "a6_1", stop_not_found=False)        
        new_lines=[]
        if epsilon_exists == 0:
            new_lines.append("                epsilon_0             Jeffreys         0.00500          0.001000          0.010000\n")
            new_lines.append("                epsilon_1             Fix              0.000000\n")
        if epsilon_exists == 1 or epsilon_exists >2:
            raise ValueError("Error: Could not find the correct number of epsilon parameters in source file. Manual debug required.")
        if theta0_exists == False:
            new_lines.append("                   theta0             Uniform          50.00000          0.000000          90.00000\n")
        if delta_exists == False:
            new_lines.append("                    delta             Uniform          20.00000          0.000000          45.00000\n")
        lines = lines[:index+1] + new_lines + lines[index+1:]
    # Find index of line with parameter 'a1_0' so that we can put the new lines just after a1_0
    index=search_index(lines, "a3_1", stop_not_found=False)           

    # Find the index of Visibility_l3 to place extra lines, if required
    Height_exists=False
    Width_exists=False
    trunc_exists=False
    for i, line in enumerate(lines[:]):
        words= line.split() # Look for only words separated by spaces or tab... This avoid True when eg. model_MS_Global_a1etaa3_HarveyLike
        if "Height" in words:
            Height_exists=True
        if "Width" in words:
            Width_exists=True
        if "trunc_c" in words:
            trunc_exists=True
    if Height_exists == False:
        lines.append("                   Height             Jeffreys          1.000000          1000.000\n")
    if Width_exists == False:
        lines.append("                    Width             Fix_Auto          1.000000\n")
    if trunc_exists == False:
        lines.append("                  trunc_c             Fix               30.0000\n")
    #
    return lines

def update_global_params_v0(lines, key, priorname=None, guess=None, prior_vals=None):
    '''
    A function that gets a model file in the form of an array of lines and
    updates the specified values as per the optional parameters.
    guess: update the initial guess (starting point of the MCMC)
    priorname: update the name of the prior (e.g., Fix)
    prior_vals: update the vector of values shown for the prior
    This program keeps the length of each lines as it is. The non-v0 version fix the length
    to a specified size.
    '''
    # Start the conditional loop
    content_start = False
    for i, line in enumerate(lines[:]):
        # No need to bother with what is before this statement...
        if line.startswith("# Controls and priors for common parameters"):
            content_start = True
        if content_start:
            if line.strip().startswith(key):
                parts = line.strip().split() # Get each field value
                req_parts = re.findall(r'\s*\S+\s*', line.strip()) # Get the size + blanks by striping only the right side ==> get the req lenghts
                # to count the number of spaces before the first character
                c=""
                for e in line.split(" "):
                    if e == "":
                        c=c + "".join(" ")
                    else: # We stop at the first non-white space
                        break
                req_parts[0]="".join(c) + parts[0]
                try:
                    req_len_priorvals=len(req_parts[3])
                except: # If the parameter is "Fix" then no priors will be there... we need to acount for that
                    req_len_priorvals=16 # Set the default expected size anyway 
                if priorname is not None:
                    parts[1] = priorname
                if guess is not None:
                    parts[2] = str(guess)
                if prior_vals is not None:
                    parts=parts[0:3]  
                    for pr in prior_vals:
                        parts.append(str(pr))
                # Ensure that we have the correct padding for all of the elements, even for the Key
                for k in range(len(parts)):
                    try:
                        parts[k]=format_string_length(parts[k], len(parts[k]), len(req_parts[k]))
                        #print(" parts[{}] = {}".format(k, parts[k]))
                    except: # This is then the case of the replacement of a Fix parameter
                        parts[k]=format_string_length(parts[k], len(parts[k]), req_len_priorvals)
                lines[i] = "".join(parts) + "\n"
    return lines

def update_global_params(lines, key, priorname=None, guess=None, prior_vals=None):
    '''
    A function that gets a model file in the form of an array of lines and
    updates the specified values as per the optional parameters.
    guess: update the initial guess (starting point of the MCMC)
    priorname: update the name of the prior (e.g., Fix)
    prior_vals: update the vector of values shown for the prior
    '''
    # Define the lengths of each section. Assumes 7 Sections: 
    #           key / priorname / guess /4 params max for prior_vals
    len_parts=[25,30,20,25,25,25,25]
    # Start the conditional loop
    content_start = False
    #for i, line in enumerate(lines):
    for i in range(len(lines)):
        line=lines[i]
        # No need to bother with what is before this statement...
        if line.startswith("# Controls and priors for common parameters"):
            content_start = True
        if content_start:
            parts = line.strip().split() # Get each field value
            #print(" key = {}  ==> line.strip().startswith(key) = {}".format(key, line.strip().startswith(key)))
            if line.strip().startswith(key):
                if priorname is not None:
                    parts[1] = priorname
                if guess is not None:
                    parts[2] = str(guess)
                if prior_vals is not None:
                    parts=parts[0:3]  
                    for pr in prior_vals:
                        parts.append(str(pr))
                print("         Modified parts = ", parts)
            # Ensure that we have the correct padding for all of the elements, even for the Key
            if line.startswith("# Controls and priors for common parameters") != True and \
                    line.strip().startswith("model_fullname") != True:
                for k in range(len(parts)):
                    parts[k]=format_string_length(parts[k], len(parts[k]), len_parts[k])
                lines[i] = "".join(parts) + "\n"
    return lines

def reformat_v1001_Kamiaka(filein, outputdir):
    '''
    This setup was used to convert the OLD Kamiaka2018 setups into the new ones
    This was used to make a 1001 setup
    '''
    file_name = os.path.basename(filein)
     # Construct the output file path
    output_file_path = os.path.join(outputdir, file_name)
    # Readjust the format of the global parameters to be compatible to either to 
    # model_MS_Global_aj_HarveyLike or model_MS_Global_ajAlm_HarveyLike
    lines=polish_file(filein, is_ajAlm_model=False)
    # Apply the rules of transformation for that version of the reformating
    lines=update_global_params(lines, "Inclination", priorname="Uniform", prior_vals=[0,90.000])
    lines=update_global_params(lines, "Visibility_l2", prior_vals=[0.53, 0.03])
    # Write the modified lines to the output file
    with open(output_file_path, 'w') as f:
        f.writelines(lines)

def reformat_Tabulated_a1inc_aj(filein, outputdir, a1_guess=None, inc_guess=None, a2_range=None, a3_range=None,a4_range=None):
    '''
    This setup was used to convert the OLD Kamiaka2018 + Legacy setups into the new ones
    This was used to make a 1111 or a 1101 setup (a1a2a3a4 fit) with Tabulated (a1,inc) prior
    filein: The file to be used as reference
    outputdir: The directory for the output file
    a1_guess and inc_guess: Allows you to update the guesses for a1 or inc. Useful to start in a non-0 Probability
        region of the Correlation MAP and thus avoid starting with a -Inf in logPosterior during the MCMC (which stops the process)
    a2_range, a3_range, a4_range: Ranges (min, max) for the uniform prior of these parameters
        If not provided, these parameters will be Fixed to 0
    '''
    is_ajAlm_model=False
    file_name = os.path.basename(filein)
     # Construct the output file path
    output_file_path = os.path.join(outputdir, file_name)
    # Readjust the format of the global parameters to be compatible to either to 
    # model_MS_Global_aj_HarveyLike or model_MS_Global_ajAlm_HarveyLike
    lines=polish_file(filein, is_ajAlm_model=is_ajAlm_model)
    # Apply the rules of transformation for that version of the reformating
    if a1_guess == None:
        lines=update_global_params(lines, "a1_0", priorname="Tabulated_2d(Inclination)", prior_vals=[0])
    else:
        lines=update_global_params(lines, "a1_0", guess=a1_guess, priorname="Tabulated_2d(Inclination)", prior_vals=[0])

    if a2_range is not None:
        lines=update_global_params(lines, "a2_0", guess=0, priorname="Uniform", prior_vals=a2_range)
    else:
        lines=update_global_params(lines, "a2_0", priorname="Fix", guess=0, prior_vals=[])
    if a3_range is not None:
        lines=update_global_params(lines, "a3_0", guess=0, priorname="Uniform", prior_vals=a3_range)
    else:
        lines=update_global_params(lines, "a3_0", priorname="Fix", guess=0, prior_vals=[])
    if a4_range is not None:
        lines=update_global_params(lines, "a4_0", guess=0, priorname="Uniform", prior_vals=a4_range)
    else:
        lines=update_global_params(lines, "a4_0", priorname="Fix", guess=0, prior_vals=[])
    if inc_guess == None:
        lines=update_global_params(lines, "Inclination", priorname="Uniform", prior_vals=[0,90.000])
    else:
        lines=update_global_params(lines, "Inclination", guess=inc_guess, priorname="Uniform", prior_vals=[0,90.000])

    # Write the modified lines to the output file
    with open(output_file_path, 'w') as f:
        f.writelines(lines)

def reformat_Tabulated_a1inc_ajAlm(filein, outputdir, a1_guess=None, inc_guess=None, a3_range=None, filter_type="gate", decompose="-1"):
    '''
    This setup was used to convert the OLD Kamiaka2018 + Legacy setups into the new ones
    This was used to make a 1211 or a 1201 setup (a1,activity,a3 fit) with Tabulated (a1,inc) prior
    filein: The file to be used as reference
    outputdir: The directory for the output file
    a1_guess and inc_guess: Allows you to update the guesses for a1 or inc. Useful to start in a non-0 Probability
        region of the Correlation MAP and thus avoid starting with a -Inf in logPosterior during the MCMC (which stops the process) 
    a3_range: Ranges (min, max) for the uniform prior of these parameters
        If not provided, these parameters will be Fixed to 0
    '''
    is_ajAlm_model=True
    file_name = os.path.basename(filein)
     # Construct the output file path
    output_file_path = os.path.join(outputdir, file_name)
    # Readjust the format of the global parameters to be compatible to either to 
    # model_MS_Global_aj_HarveyLike or model_MS_Global_ajAlm_HarveyLike
    lines=polish_file(filein, is_ajAlm_model=is_ajAlm_model)
    # Apply the rules of transformation for that version of the reformating
    if a1_guess == None:
        lines=update_global_params(lines, "a1_0", priorname="Tabulated_2d(Inclination)", prior_vals=[0])
    else:
        lines=update_global_params(lines, "a1_0", guess=a1_guess, priorname="Tabulated_2d(Inclination)", prior_vals=[0])
    lines=update_global_params(lines, "filter_type", guess=None, priorname=filter_type, prior_vals=None)
    lines=update_global_params(lines, "decompose_Alm", guess=decompose, priorname="Fix", prior_vals=None)
    if a3_range is not None:
        lines=update_global_params(lines, "a3_0", guess=0, priorname="Uniform", prior_vals=a3_range)
    else:
        lines=update_global_params(lines, "a3_0", priorname="Fix", guess=0, prior_vals=[])
    if inc_guess == None:
        lines=update_global_params(lines, "Inclination", priorname="Uniform", prior_vals=[0,90.000])
    else:
        lines=update_global_params(lines, "Inclination", guess=inc_guess, priorname="Uniform", prior_vals=[0,90.000])
    # Write the modified lines to the output file
    with open(output_file_path, 'w') as f:
        f.writelines(lines)

def process_stars(dir_in, dir_priors, dir_out_root, a2a3a4_table):
    '''
        Use tables of a2,a3,a4, reference model files and priors tables to prepare the configurations
    '''
    dir_out_1111= dir_out_root + "/1111/"
    dir_out_1101= dir_out_root + "/1101/"
    dir_out_1011= dir_out_root + "/1011/"
    dir_out_1211= dir_out_root + "/1211/"
    dir_out_1201= dir_out_root + "/1201/"
    for data in a2a3a4_table:
        # Preparing required information
        identifier=data[0]
        filein=find_modelfile_byID(dir_in, identifier, prefix="kplr", suffix="", extension=".model")
        filepriors=find_modelfile_byID(dir_priors, identifier, prefix="kplr", suffix="", extension=".priors")

        a1_guess, inc_guess, Proba=find_max_prior(dir_priors + filepriors)
        print("      a1_guess = ", a1_guess)
        print("      inc_guess = ", inc_guess)
        print(" Processing: ", dir_in+filein)
        print("     Case 1111:")
        reformat_Tabulated_a1inc_aj(dir_in+filein, dir_out_1111, a1_guess=a1_guess, inc_guess=inc_guess, a2_range=data[2:4], 
                                    a3_range=data[4:6], a4_range=data[6:])
        print("     Case 1101:")
        reformat_Tabulated_a1inc_aj(dir_in+filein, dir_out_1101, a1_guess=a1_guess, inc_guess=inc_guess,  a2_range=data[2:4], 
                                    a3_range=None, a4_range=data[6:])
        print("     Case 1011:")
        reformat_Tabulated_a1inc_aj(dir_in+filein, dir_out_1011, a1_guess=a1_guess, inc_guess=inc_guess,  a2_range=None, 
                                    a3_range=data[4:6], a4_range=None)
        # ----
        # ----
        # ----
        print("     Case 1211:")
        print("  ==> Gate -> decompose = -1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1211 + "/Gate/decompose_-1/", 
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=data[4:6], 
                                       filter_type="gate", decompose="-1")
        print("  ==> Gate -> decompose = 1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1211 + "/Gate/decompose_1/", 
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=data[4:6], 
                                       filter_type="gate", decompose="1")
        print("  ==> Gate -> decompose = 2")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1211 + "/Gate/decompose_2/", 
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=data[4:6], 
                                       filter_type="gate", decompose="2")
        print("  ==> Triangle -> decompose = -1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1211 + "/Triangle/decompose_-1/", 
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=data[4:6], 
                                       filter_type="triangle", decompose="-1")
        print("  ==> Triangle -> decompose = 1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1211 + "/Triangle/decompose_1/", 
                                        a1_guess=a1_guess, inc_guess=inc_guess, a3_range=data[4:6], 
                                       filter_type="triangle", decompose="1")
        print("  ==> Triangle -> decompose = 2")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1211 + "/Triangle/decompose_2/", 
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=data[4:6], 
                                       filter_type="triangle", decompose="2")
        # ---- 
        # ----
        # ----
        print("     Case 1201:")
        print("  ==> Gate -> decompose = -1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1201 + "/Gate/decompose_-1/",
                                       a1_guess=a1_guess, inc_guess=inc_guess,  a3_range=None, 
                                       filter_type="gate", decompose="-1")
        print("  ==> Gate -> decompose = 1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1201 + "/Gate/decompose_1/",
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=None, 
                                       filter_type="gate", decompose="1")
        print("  ==> Gate -> decompose = 2")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1201 + "/Gate/decompose_2/",
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=None, 
                                       filter_type="gate", decompose="2")
        print("  ==> Triangle -> decompose = -1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1201 + "/Triangle/decompose_-1/",
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=None, 
                                       filter_type="triangle", decompose="-1")
        print("  ==> Triangle -> decompose = 1")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1201 + "/Triangle/decompose_1/",
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=None, 
                                       filter_type="triangle", decompose="1")
        print("  ==> Triangle -> decompose = 2")
        reformat_Tabulated_a1inc_ajAlm(dir_in+filein, dir_out_1201 + "/Triangle/decompose_2/",
                                       a1_guess=a1_guess, inc_guess=inc_guess, a3_range=None, 
                                       filter_type="triangle", decompose="2")      
    print(" ------ model files were created in directories with this structure -----")
    margin=" "*5
    print(" root directory: ", dir_out_root)
    print(" subdirectory structure:")
    print(margin, " - 1111")
    print(margin, " - 1011")    
    print(margin, " - 1101")
    print(margin, "                     --> decompose_-1")
    print(margin, "                    |")
    print(margin, " - 1201 ------> Gate --> decompose_1")
    print(margin, "       |            |")
    print(margin, "       |             --> decompose_2")
    print(margin, "       |")
    print(margin, "       |             --> decompose_-1")
    print(margin, "       |            |")
    print(margin, "        --> Triangle --> decompose_1")
    print(margin, "                    |")
    print(margin, "                     --> decompose_2")
    print("\n")
    print(margin, "                     --> decompose_-1")
    print(margin, "                    |")
    print(margin, " - 1211 ------> Gate --> decompose_1")
    print(margin, "       |            |")
    print(margin, "       |             --> decompose_2")
    print(margin, "       |")
    print(margin, "       |             --> decompose_-1")
    print(margin, "       |            |")
    print(margin, "        --> Triangle --> decompose_1")
    print(margin, "                    |")
    print(margin, "                     --> decompose_2")


def do_setups_legacy():
    '''
        Main code to derive setups for LEGACY stars and for the analysed stars and based on the initial 
        run 1001
    '''
    # -------------------------
    # ------ LEGACY STARS -----
    # -------------------------
    # Table of a2,a3,a4 inputs for the Legacy stars
    file_table_legacy="../../data/REFERENCE_RUN_1001/Legacy_stars/a2a3a4_priortable_Legacy.txt"
    _,a2a3a4_table=read_a2a3a4_prior_table(file_table_legacy)

    # Location of the 2D priors tables. Used to update the (a1,inc) guess to be sure that
    # we start in a non-Inf Probability value (thus starting with a better guess)
    dir_priors="../../data/TRANSFORMED_RUN/Inputs/Legacy/priors/"

    # Reference Inputs
    dir_in="../../data/REFERENCE_RUN_1001/Legacy_stars/inputs_1001_science2018/1001/"
    
    # Output Directories
    dir_out_root= "../../data/TRANSFORMED_RUN/Inputs/Legacy/"
    process_stars(dir_in, dir_priors, dir_out_root, a2a3a4_table)


def do_setups_kamiaka():
    '''
        Main code to derive setups for Kamiaka stars for the analysed stars and based on the initial 
        run 1001
    '''
    # -----------------------------------------
    # ------ KAMIAKA STARS WITHOUT LEGACY -----
    # -----------------------------------------
    # Table of a2,a3,a4 inputs for the Legacy stars
    file_table_legacy="../../data/REFERENCE_RUN_1001/Kamiaka_stars/a2a3a4_priortable_Kamiaka.txt"
    _,a2a3a4_table=read_a2a3a4_prior_table(file_table_legacy)

    # Location of the 2D priors tables. Used to update the (a1,inc) guess to be sure that
    # we start in a non-Inf Probability value (thus starting with a better guess)
    dir_priors="../../data/TRANSFORMED_RUN/Inputs/Kamiaka2018/26Feb2024/priors/"

    # Reference Inputs. Those are 1001, reprocessed Kamiaka stars
    #dir_in="../../data/REFERENCE_RUN_1001/Kamiaka_stars/14Aug2023/inputs/"
    dir_in="../../data/REFERENCE_RUN_1001/Kamiaka_stars/26Feb2024/inputs/1001/"

    # Output Directories for the model and data files that will be created
    dir_out_root= "../../data/TRANSFORMED_RUN/Inputs/Kamiaka2018/26Feb2024/"

    process_stars(dir_in, dir_priors, dir_out_root, a2a3a4_table)

def redo_1001_legacy():
    '''
        This is to redo a setup for 1001 based on the Legacy 1001
        This new setup would contain the 2D prior on a1 and inclination. This to allow 
        a comparison with other runs that do include this priors and also to deal with
        a potential issue in the normalisation constant that I noted for the 2D prior. 
        
    '''
    # -------------------------
    # ------ LEGACY STARS -----
    # -------------------------
    # Table of a2,a3,a4 inputs for the Legacy stars
    file_table_legacy="../../data/REFERENCE_RUN_1001/Legacy_stars/a2a3a4_priortable_Legacy.txt"
    _,a2a3a4_table=read_a2a3a4_prior_table(file_table_legacy)

    # Location of the 2D priors tables. Used to update the (a1,inc) guess to be sure that
    # we start in a non-Inf Probability value (thus starting with a better guess)
    dir_priors="../../data/TRANSFORMED_RUN/Inputs/Legacy/priors/"

    # Reference Inputs
    dir_in="../../data/REFERENCE_RUN_1001/Legacy_stars/inputs_1001_science2018/1001/"
    
    # Output Directories
    dir_out_root= "../../data/TRANSFORMED_RUN/Inputs/Legacy/"
    dir_out_1001= dir_out_root + "/1001/"

    print("     Case 1001 (redo):")
    for data in a2a3a4_table:
        # Preparing required information
        identifier=data[0]
        filein=find_modelfile_byID(dir_in, identifier, prefix="kplr", suffix="", extension=".model")
        filepriors=find_modelfile_byID(dir_priors, identifier, prefix="kplr", suffix="", extension=".priors")

        a1_guess, inc_guess, Proba=find_max_prior(dir_priors + filepriors)
        print("      a1_guess = ", a1_guess)
        print("      inc_guess = ", inc_guess)
        reformat_Tabulated_a1inc_aj(dir_in+filein, dir_out_1001, a1_guess=a1_guess, inc_guess=inc_guess,  a2_range=None, 
                                    a3_range=None, a4_range=None)

def redo_1001_kamiaka():
    '''
        This is to redo a setup for 1001 based on the Kamiaka 1001
        This new setup would contain the 2D prior on a1 and inclination. This to allow 
        a comparison with other runs that do include this priors and also to deal with
        a potential issue in the normalisation constant that I noted for the 2D prior. 
        
    '''
    # -----------------------------------------
    # ------ KAMIAKA STARS WITHOUT LEGACY -----
    # -----------------------------------------
    # Table of a2,a3,a4 inputs for the Legacy stars
    file_table_legacy="../../data/REFERENCE_RUN_1001/Kamiaka_stars/a2a3a4_priortable_Kamiaka.txt"
    _,a2a3a4_table=read_a2a3a4_prior_table(file_table_legacy)

    # Location of the 2D priors tables. Used to update the (a1,inc) guess to be sure that
    # we start in a non-Inf Probability value (thus starting with a better guess)
    dir_priors="../../data/TRANSFORMED_RUN/Inputs/Kamiaka2018/26Feb2024/priors/"

    # Reference Inputs. Those are 1001, reprocessed Kamiaka stars
    #dir_in="../../data/REFERENCE_RUN_1001/Kamiaka_stars/14Aug2023/inputs/"
    dir_in="../../data/REFERENCE_RUN_1001/Kamiaka_stars/26Feb2024/inputs/1001/"

    # Output Directories
    dir_out_root= "../../data/TRANSFORMED_RUN/Inputs/Kamiaka2018/26Feb2024/"
    dir_out_1001= dir_out_root + "/1001/"

    print("     Case 1001 (redo):")
    for data in a2a3a4_table:
        # Preparing required information
        identifier=data[0]
        filein=find_modelfile_byID(dir_in, identifier, prefix="kplr", suffix="", extension=".model")
        filepriors=find_modelfile_byID(dir_priors, identifier, prefix="kplr", suffix="", extension=".priors")

        a1_guess, inc_guess, Proba=find_max_prior(dir_priors + filepriors)
        print("      a1_guess = ", a1_guess)
        print("      inc_guess = ", inc_guess)
        reformat_Tabulated_a1inc_aj(dir_in+filein, dir_out_1001, a1_guess=a1_guess, inc_guess=inc_guess,  a2_range=None, 
                                    a3_range=None, a4_range=None)
