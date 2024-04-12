'''
    Programs in charge of:
        - getting the Evidence for multiple models
        - Extract the most likely among user-specified list of models 
        - Make an image of a matrix reprensenting all of the evidences for given models
'''
import os
import glob
import numpy as np
from pythonlibs.read_outputs_tamcmc import read_global_likelihood
from pythonlibs.read_outputs_tamcmc import getevidence
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from termcolor import colored

def compute_evidence_by_models(data_dir, process_name, phase="A", prefix="kplr", cpp_path="cpp_prg/"):
    '''
        This function (re)compute the evidences using the getevidence function available since 
        version 1.86.75 of TAMCMC
        It is calling read_outputs_tamcmc.py::getevidence() with a preset of parameters
        Args:
            data_dir (str): Root directory to search for subdirectories.
            phase (str): Phase identifier (default is "A").
            prefix (str): Optional prefix to filter subdirectories (default is None).
            cpp_path (str): Optional location directory for the compiled C++ program 'getevidence'
    '''
    extension_name='_evidence.txt' # We overwrite the existing file
    evidence, err_evidence=getevidence(data_dir, process_name, phase=phase, interpolation_factor=1000, 
                first_index=0, last_index=-1, period=7, bootstrap_method="block",
                bootstrap_Nsamples=2000, boostrap_blocksize_percent=1,
                erase_tmp=True, cpp_path=cpp_path, extension_name=extension_name, outdir=None)
    return evidence, err_evidence

def test_getevidence():
    data_dir="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Outputs/Tcoef_1/Legacy/1001/"
    process_name="kplr008760414_kasoc-psd_slc_v1"
    phase="A"
    cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.86.75_AMD/"
    evidence, err_evidence=compute_evidence_by_models(data_dir, process_name, phase=phase, cpp_path=cpp_path)
    print("evidence : ", evidence, flush=True)
    print("err_evidence:", err_evidence, flush=True)

def get_evidence_by_models(data_dir, phase="A", prefix="kplr", recompute=[False, False]):
    '''
        This function read all of the evidence stored in the outputs directories of the MCMC processes
        Given a root directory data_dir, go into all of the analysed stars and 
        retrieve the evidence for each of them. 
        Output results are stored in a numpy arrays (for the evidence and its error)
        and one list for the star ID
        
        Args:
            data_dir (str): Root directory to search for subdirectories.
            phase (str): Phase identifier (default is "A").
            prefix (str): Optional prefix to filter subdirectories (default is None).
            recompute (list): Optional list that defines if the evidence is recomputed using getevidence (available since TAMCMC 1.86.75)
                              Minimal content is recompute=[False, False] for no recomputation and bootstraping OR [False, True] for no recomputation and results from bootstraping
                              Otherwise, recompute=[True, True, cpp_path]. Note that the bootstrap must be on
                              bootstrap_on (bool): Optional argument specifying if the file to was computed using bootstrap_on=True or False
                              If you read directly outputs of the TAMCMC without using getevidence, you should set bootstrap_on=False
                              bootstraping computation being a bit slow, it is not performed iteratively when writting data from the TAMCMC
        Returns:
            tuple: Numpy arrays for evidence and its error, and a list for star IDs.
    '''
    subdirectories = [name for name in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, name)) and (prefix is None or name.startswith(prefix))]
    
    evidence_list = []
    err_evidence_list = []
    star_id_list = []
    extension_name="evidence.txt"
    for subdirectory in subdirectories:
        file_in = os.path.join(data_dir, subdirectory, "diags", f"{subdirectory}_{phase}_"+extension_name)
        if recompute[0] == False:
            # Check if the file exists
            if not os.path.exists(file_in):
                print(f" Warning: File {file_in} does not exist. Attempting to find a file with a more relaxed syntax...", flush=True) 
                # Try to find a file with the relaxed syntax
                relaxed_file_pattern = os.path.join(data_dir, subdirectory, "diags", f"{subdirectory}_*_{phase}_evidence.txt")
                matching_files = glob.glob(relaxed_file_pattern)   
                # If a matching file is found, use it as the new file_in
                if matching_files:
                    file_in = matching_files[0]
                    print(f"      Alternative file found... Using file: {file_in}", flush=True)
                else:
                    print("Error finding evidence file. Try to recompute it", flush=True)
                    raise FileNotFoundError(f"No file found with the syntax: {relaxed_file_pattern}")
            evidence, err_evidence = read_global_likelihood(file_in, evidence_only=True, bootstrap_on=recompute[1])
        else:
            evidence, err_evidence=compute_evidence_by_models(data_dir, subdirectory, phase="A", prefix="kplr", cpp_path=recompute[2])
        evidence_list.append(evidence)
        err_evidence_list.append(err_evidence)
        star_id_list.append(subdirectory)
    return np.asarray(evidence_list), np.asarray(err_evidence_list), star_id_list

def translate_modeltype2dir(root_dir, modeltype, check_dir_exist=True, no_dir_is_fatal=True):
    '''
        Given a root directory that contain all of the models computed following a tree structure
        reconstruct the location of all the stellar data of a given model.
        for example : type = [1101] ,
        Should give, dir_out=root_dir + "1101/"
        and  type = [1201, "Triangle", "decompose-1"]
        Should give, dir_out=root_dir + "1201/Triangle/decompose-1/"
             " - 1111"
             " - 1101"
             "                     --> decompose_-1"
             "                    |"
             " - 1201 ------> Gate --> decompose_1"
             "       |            |"
             "       |             --> decompose_2"
             "       |")
             "       |             --> decompose_-1"
             "       |            |"
             "        --> Triangle --> decompose_1"
             "                    |"
             "                     --> decompose_2"
             
             "                     --> decompose_-1"
             "                    |")
             " - 1211 ------> Gate --> decompose_1"
             "       |            |")
             "       |             --> decompose_2"
             "       |")
             "       |             --> decompose_-1"
             "       |            |")
             "        --> Triangle --> decompose_1"
             "                    |")
             "                     --> decompose_2"
            '''
    dir_out=None
    if len(modeltype) == 1:
        dir_out = root_dir + "/" + str(modeltype[0])
    elif len(modeltype) == 3:
        dir_out = root_dir + "/" + str(modeltype[0]) + "/" + str(modeltype[1]) + "/" + str(modeltype[2])
    else:
        raise ValueError("Invalid modeltype format")
    if check_dir_exist==True:
        if not os.path.exists(dir_out):
            if no_dir_is_fatal == True:
                raise ValueError("Directory {} does not exist".format(dir_out))
            else:
                return None    
    return dir_out

def translate_dir2modeltype(root_dir, full_path):
    '''
    Given a root directory and a full path, extract the modeltype from the full path.
    '''
    path_parts = full_path.replace(root_dir, "").split("/")
    modeltype = [int(part) if part.isdigit() else part for part in path_parts if part]
    allowed_models = ["1001", "1011","1101", "1111", "1201", "1211"]
    if str(modeltype) not in allowed_models:
        raise ValueError("Not supported models. You must set them to 1001, 1101, 1111, 1201 or 1211")
    return modeltype

def compute_Odds_ratios(root_dir, phase="A", prefix="kplr", 
                            allowed_modeltype=[1001, 1101,1111,1201,1211],
                            allowed_modelsubtype1=["Gate", "Triangle"],
                            allowed_modelsubtype2=["decompose_-1", "decompose_1", "decompose_2"],
                            return_evidences=False, recompute_evidence=[False], reference_dir="1001"
                        ):
    # Define allowed types of models by generating combinations of all these type and subtypes
    modeltypes = []
    if any("1001" in str(model) for model in allowed_modeltype):
        modeltypes.append([1001])
    if any("1101" in str(model) for model in allowed_modeltype):
        modeltypes.append([1101])
    if any("1011" in str(model) for model in allowed_modeltype):
        modeltypes.append([1011])
    if any("1111" in str(model) for model in allowed_modeltype):
        modeltypes.append([1111])
    for modeltype in allowed_modeltype:
        for subtype1 in allowed_modelsubtype1:
            for subtype2 in allowed_modelsubtype2:
                if modeltype not in [1001, 1011,1101, 1111]:
                    modeltypes.append([modeltype, subtype1, subtype2])
    
    # Get the number of stars by reading the reference directory
    # WE WILL ASSUME THAT ALL OTHER DIRECTORIES HAVE THE SAME NUMBER OF STARS
    model_ref_dir=os.path.join(root_dir, reference_dir)
    subdirectories = [name for name in os.listdir(model_ref_dir) if os.path.isdir(os.path.join(model_ref_dir, name)) and (prefix is None or name.startswith(prefix))]
    Nstars=len(subdirectories)
    # Generate all of the possible model combinations 
    oddsratios_all=np.zeros((Nstars, len(modeltypes), len(modeltypes)))
    oddsratios_err_all=np.zeros((Nstars,len(modeltypes), len(modeltypes)))

    for i in range(len(modeltypes)):
        for j in range(len(modeltypes)):
            dir_modelA=translate_modeltype2dir(root_dir, modeltypes[i], check_dir_exist=True)
            model_evidencesA, model_err_evidencesA, star_id_listA=get_evidence_by_models(dir_modelA, phase=phase, recompute=recompute_evidence)
            dir_modelB=translate_modeltype2dir(root_dir, modeltypes[j], check_dir_exist=True)
            model_evidencesB, model_err_evidencesB, star_id_listB=get_evidence_by_models(dir_modelB, phase=phase, recompute=recompute_evidence)
            if star_id_listA != star_id_listB:
                print(" Nstars =", Nstars, flush=True)
                print(" len(star_id_listA) =", len(star_id_listA), flush=True)
                print(" len(star_id_listB) =", len(star_id_listB), flush=True)
                for i in range(len(star_id_listA)):
                    print(" {}   {} ".format(star_id_listA[i], star_id_listB[i]) , flush=True)
                raise ValueError("Error: Lists of stars found to be different while computing the evidences.")
            odds=np.exp(model_evidencesA-model_evidencesB) # >1 means P(A) > P(B) or : in i > in j
            oddsratios_all[:,i,j]=odds[:]
            err_odds=(model_err_evidencesA**2 + model_err_evidencesB**2)
            err_odds=np.sqrt(err_odds) * np.exp(model_evidencesA-model_evidencesB)
            oddsratios_err_all[:,i,j]=err_odds
    if return_evidences == False:
        return modeltypes, oddsratios_all, oddsratios_err_all, star_id_listA
    else:
        return modeltypes, oddsratios_all, oddsratios_err_all, star_id_listA, [], []
    
def make_short_names(modeltype):
    passed=False
    if str(modeltype[0]) == "1001":
        passed=True
        return "1001"
    if str(modeltype[0]) == "1011":
        passed=True
        return "1011"
    if str(modeltype[0]) == "1101":
        passed=True
        return "1101"
    if str(modeltype[0]) == "1111":
        passed=True
        return "1111"
    if str(modeltype[0]) == "1201" or str(modeltype[0]) == "1211":
        passed=True
        end=modeltype[2].split("_")
        return str(modeltype[0]) + modeltype[1][0] + end[-1]
    if passed == False:
        raise ValueError("Unidentified modeltype: ", modeltype)

def show_probabilities(modeltypes, oddsratios, oddsratios_err, fileout, starID):
    '''
    Receives a square matrix of Probability and make a plot 
    showing the associated probability of each element.
    Also return a ranked probability in a list, with element 0 
    being the highest probability one
    '''
    # Create color map for background colors
    cmap = plt.cm.get_cmap('RdYlGn')
    norm = LogNorm(vmin=0.001, vmax=oddsratios.max())

    np.seterr(divide='raise')
    # Create the plot
    fig, ax = plt.subplots(1, figsize=(12, 6), num=1, clear=True)
    im=ax.imshow(oddsratios, cmap=cmap, norm=norm)

    probabilities=np.zeros(len(starID))
    probabilities_min=np.zeros(len(starID))
    probabilities_max=np.zeros(len(starID))
    
    # Add probability values and error bounds to each element
    for i in range(len(modeltypes)):
        try:
            pi = 100 / (np.sum(1. / oddsratios[i,:]))
            plist=[100 / (np.sum(1. / (oddsratios[i,:]-oddsratios_err[i,:]))), 100 / (np.sum(1. / (oddsratios[i,:]+oddsratios_err[i,:])))]
        except FloatingPointError:
            print("Pi =", pi, flush=True)
            print(" Oddsratios = ", oddsratios[i,:], flush=True)
            exit()
        probabilities[i]=pi
        probabilities_min[i]=np.min(plist)
        probabilities_max[i]=np.max(plist)
        for j in range(len(modeltypes)):
            text = ax.text(j, i, f'{oddsratios[i, j]:.1e}' if oddsratios[i, j] > 100 else f'{oddsratios[i, j]:.1f}', ha='center', va='center', fontsize=4)

    # Add the last column for probabilities
    for i in range(len(modeltypes)):
        text = ax.text(len(modeltypes), i, f'{probabilities[i]:.1f}', ha='center', va='center', fontsize=4, color='black')

    tick_names = []
    for m in modeltypes:
        tick_names.append(make_short_names(m))

    # Set axis labels and title
    ax.set_xticks(np.arange(len(tick_names) + 2))
    ax.set_yticks(np.arange(len(tick_names)))
    ax.set_xticklabels(tick_names + ['Probabilities',''])
    ax.set_yticklabels(tick_names)
    ax.set_title(starID)

    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Show the plot
    fig.savefig(fileout, dpi=400)
    plt.close()
    return modeltypes, probabilities, probabilities_min, probabilities_max


def compile_results_comparison(starID, modeltypes, oddsratios, oddsratios_err, dir_out, 
                               file_summary="odds_summary", append=False, 
                               color_thld=[75,95], withcolors=True):
    '''
        Perform the comparison of probability and save it as an image as well as a table
        starID: List of stars Identifications
        modeltypes: List of type of model analysed
        oddsratios: Array of oddratios
        oddsratios: Array of errors on the oddratios
        dir_out: output directory
        file_summary: Name of the output file
        append: If False, create a new file. Otherwise append it
        withcolors: If True, use colors to highlight some solutions. Deactivate when writting on files
        color_thld: Defines the thresholds from which we highlight either in green, yellow or red the Probability
    '''
    # Iterate over each starID
    i=0
    for s in range(len(starID)):
        if append == True:
            mode="a"
        else:
            mode="w"
        # Print the starID to the console
        print("Processing star: ", starID[s], flush=True)
        # Create the file path for the output image
        fileout_jpg = os.path.join(dir_out , starID[s] + ".jpg")
        # Get the modeltypes_out and probabilities
        modeltypes_out, probabilities, probabilities_min, probabilities_max = show_probabilities(modeltypes, oddsratios[s,:,:], oddsratios_err[s,:,:], fileout_jpg, starID[s])
        # Initialize the strings for modeltypes and probabilities
        str_model = "("
        string_med = ""
        string_m1s = ""
        string_p1s = ""
        # Iterate over each modeltype and probability
        for j in range(len(modeltypes_out)):
            # Determine the color based on the probability value
            if withcolors == True:
                if probabilities[j] < color_thld[0]:
                    col = "red"
                elif probabilities[j] >= color_thld[0] and probabilities[j] <= color_thld[1]:
                    col = "yellow"
                else:
                    col = "green"
                # Append the formatted modeltype and probability to the strings
                str_model += colored("{}".format(str(modeltypes_out[j])), col)
                string_med += colored("{:<26.1f}".format(probabilities[j]), col)
                string_m1s += colored("{:<26.1f}".format(probabilities_min[j]), col)
                string_p1s += colored("{:<26.1f}".format(probabilities_max[j]), col)
            else:
                # Append the formatted modeltype and probability to the strings
                str_model += "{}".format(str(modeltypes_out[j]))
                string_med += "{:<26.1f}".format(probabilities[j])
                string_m1s += "{:<26.1f}".format(probabilities_min[j])
                string_p1s += "{:<26.1f}".format(probabilities_max[j])
        # Print the modeltypes and probabilities to the console
        print(str_model, ") :", string_med, flush=True)
        print('--', flush=True)
        # Write the starID, modeltype, and probability to the file
        file_out=file_summary + "_median.txt"
        with open(os.path.join(dir_out, file_out), mode) as file:
            if i == 0 and append == False:
                str_title="{0:<50s}".format("starID")
                for j in range(len(modeltypes_out)):
                    s_m=""
                    for m in modeltypes_out[j]:
                        s_m = s_m + str(m)
                    str_title=str_title + "{:<26s}".format(s_m)
                file.write(str_title +"\n")
            file.write("{:<50s} {}\n".format(starID[s],  string_med))
        file_out=file_summary + "_m1s.txt"
        with open(os.path.join(dir_out, file_out), mode) as file:
            if i == 0 and append == False:
                str_title="{0:<50s}".format("starID")
                for j in range(len(modeltypes_out)):
                    s_m=""
                    for m in modeltypes_out[j]:
                        s_m = s_m + str(m)
                    str_title=str_title + "{0:<26s}".format(s_m)
                file.write(str_title +"\n")
            file.write("{:<50s} {}\n".format(starID[s], string_m1s))
        file_out=file_summary + "_p1s.txt"
        with open(os.path.join(dir_out, file_out), mode) as file:
            if i == 0 and append == False:
                str_title="{0:<50s}".format("starID")
                for j in range(len(modeltypes_out)):
                    s_m=""
                    for m in modeltypes_out[j]:
                        s_m = s_m + str(m)
                    str_title=str_title + "{0:<26s}".format(s_m)
                file.write(str_title +"\n")
            file.write("{:<50s} {}\n".format(starID[s], string_p1s))
        append=True
        i=i+1

