from pythonlibs.a2CF_a2AR import getdir_idlist
from pythonlibs.a2CF_a2AR import do_a2CF_a2AR
from pythonlibs.compute_stats import *
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import re
import copy
from termcolor import colored
from pythonlibs.whisker import Draw_BoxAndWiskers_vertical
from pythonlibs.read_generic_files import read_matrix_tab

def extract_ID_from_processname(ProcessName):
    '''
        Process names often contain the star ID in it in the form of a numerical sequence
        This function extract that numerical sequence and attribute it to the star ID
    '''
    match = re.search(r'\d+', ProcessName)
    if match:
        star_id = match.group()
    else:
        # If no match is found, set the star ID as None
        star_id = None    
    return star_id

def format_ID(ID, Ndigits=9):
	'''
		Small function that ensure that the ID number is on Ndigits digits
		If this is not the case, we fill with 0 before the existing ID
		With Kepler data, Ndigits=9 is fine. For TESS, Ndigits=10 shoud be OK (need check)
	'''
	NID=len(ID)
	delta=8-NID+1 # We want to format using 8 digits
	New_ID=ID
	for d in range(delta):
		New_ID='0' + New_ID
	return New_ID

def get_externaldata(external_file):
    '''
    This function regroups the information about Mass, Radius, Teff
    The file will be returned as a list of lists due to the fact that such a file is expectected to have 
    a mixture of string and float
    '''
    #print('Reading the external data file with Teff, M, R, Age...')
    data, header, label, unit=read_matrix_tab(external_file, separator=",", first_data_is_label=True, datatype="string")
    NonSeismic={"label": label[1:], # We need to remove the starID as it is put outside now
            "header": header,
            "unit":unit[1:]
            }
    # Filtering out the starID from the rest to ease accessibility to this information
    starID=[]
    data_out=[]
    for d in data:
        if d[0] != "":
            starID.append(d[0]) # pick the KIC
            data_out.append(d[1:])
    NonSeismic["starID"]=starID
    NonSeismic["data"]=data_out
    return NonSeismic

def get_productsdata(products_root_dir, core_odds_files="Proba_summary", core_rot_files="rotation_stats",
                     odds_keys=["m1s", "median", "p1s"],
                     confidence_rot=[2.25, 16, 50, 84, 97.75]):
    '''
    Function that read synthesis files produced by do_main_summary.py
    It is assumed that the products_root_dir will have the following:
        - A "odds_ratio" subdirectory with all of the odds ratio
        - Subdirectories with names matching the entry in the data_sources list
        - A series of rotation_stats_*.summary files with the statistical summary for the rotation
    Note that this function only requires data within "odds_ratio" and the *.summary files as this contain everything.
    products_root_dir: Directory where all of the files of the postprocessing database is available
    core_odds_files: The core name for odds ratio output files. It will be appended by the strings given in the tag
                     Typically, "_m1s.txt" (minus 1 sigma confidence interval),
                     "_median.txt", "p1s.txt" (plus 1 sigma confidence interval). This will be the targeted files
    core_rot_files: The core name of rotation output files. It will be appended by the confidence_rot elements to
                    form the targeted rotation summary files.
    odds_keys: A list specifying how the odds ratio files were tagged. Used to target data files and to generate the dictionary
    confidence_rot : A list of all of the confidence intervals that were used to compute the rotation summary files
    Returns: Two dictionaries. One containing the odds ratio. The other containing the rotations.
    '''
    ProductsOdds = {"starID": {}, "Probabilities": {}, "header": None, "label": None, "Ncols":0., "Nrows":0.}  # 
    ProductsRot = {"starID":{}, "data": {}, "header": None, "label": None, "unit": None, "Ncols":0., "Nrows":0}  # 
    
    odds_ratio_dir=os.path.join(products_root_dir, "odds_ratio")
    # Read odds ratio summary files
    for i in range(len(odds_keys)):
        ofile=os.path.join(odds_ratio_dir, core_odds_files + "_"+ odds_keys[i]+".txt")
        if not os.path.exists(ofile):
            err_msg=colored("Error: Could not found the odds ratio summary file: {}".format(ofile), "red")
            raise ValueError(err_msg)
        else:
            data_odds, header, label,unit=read_matrix_tab(ofile, separator=None, first_data_is_label=True, second_data_is_unit= False, datatype="string")
            # Parse in an intelligible way the data
            ID=[]
            Pr=np.zeros((len(data_odds), len(label)-1), dtype=float)
            ProductsOdds["Nrows"]=len(data_odds)
            ProductsOdds["Ncols"]=len(label)-1
            k=0
            for d in data_odds:
                ID.append(d[0])
                Pr[k, :]=np.asarray(d[1:], dtype=float)
                k=k+1
            ProductsOdds["Probabilities"][odds_keys[i]]=Pr
            if i == 0:
                ProductsOdds["starID"]=ID
                ProductsOdds["header"]=header
                ProductsOdds["label"]=label[1:] # We remove the StarID from the list of labels
    # Read rotation summary files
    i=0
    for c in confidence_rot:
        rfile=os.path.join(products_root_dir, core_rot_files + "_{0:4.2f}.summary".format(c))
        if not os.path.exists(rfile):
            err_msg=colored("Error: Could not found the odds ratio summary file: {}".format(rfile), "red")
            raise ValueError(err_msg)
        else:
            data_rot, header, label,unit=read_matrix_tab(rfile, separator=None, first_data_is_label=True, second_data_is_unit=True, datatype="string")
            # Parse in an intelligible way the data
            ID=[]
            ModelCode=[]
            data=np.zeros((len(data_rot), len(label)-2))
            ProductsRot["Nrows"]=len(data_rot)
            ProductsRot["Ncols"]=len(label)-2
            k=0
            for d in data_rot:
                ID.append(d[0])
                ModelCode.append(d[1])
                data[k, :]=np.asarray(d[2:], dtype=float)
                k=k+1
            ProductsRot["data"][str(c)]=data
            if i == 0:
                ProductsRot["starID"]=ID
                ProductsRot["model"]=ModelCode
                ProductsRot["header"]=header
                ProductsRot["label"]=label[2:]  # We remove the StarID from the list of labels
                ProductsRot["unit"]=unit
        i=i+1
    # Perform a consistency check to verify that each 'starID' in ProductsOdds has a 'starID' counterpart in ProductsRot
    status=check_consistency(ProductsOdds, ProductsRot)
    if status == True:
        return ProductsOdds, ProductsRot
    else:
        err_msg=colored("Error in get_productsdata(): We could not find a complete match between ProductsOdds and ProductsRot\n Check visually the consistency of the data in input files!", "red")
        raise ValueError(err_msg)
# // To put in a proper test function
#products_root_dir="/Users/obenomar/Work/dev/ajAlm/test_files/statistical_summary/test/"
#ProductsOdds, ProductsRot=get_productsdata(products_root_dir)
#print("ProductsRot", ProductsRot["data"]["50"])

def query_params_from_ProductsRot(ProductsRot, ID, ModelCode, confidence, ParamName=None, ResolveID=True):
    '''
        A function that recover a specific array of parameters from ProductsRot, given the ID, the ModelCode and significance of the model
        Optionally, it can also retrieve only a single parameter using its name
        ProductsRot: The dictionary containing all of the data regarding rotation
        ID: The model identifier, either in full (in the form kplrXXXX) or in part (eg. only the KIC number)
        ModelCode: The model code associated to the model. Can be either 1001, 1101, 1111, 1201Gatedecompose-1 etc...
        confidence : The confidence level of the queried data (eg. 50 for the median and 16 for 1s lower bound)
        ParamName (Optional): the parameter name that is specifically queried (eg. "a2"). If not specified, will return the whole array of parameters
        ResolveID (Optional): If set to True (default) will keep only the ID number and remove any other parts from the ID and from within the table 
    '''
    data=[]
    if ResolveID == True:
        ID=extract_ID_from_processname(ID)  # Ensure that ID is only numbers  
        ID=format_ID(ID, Ndigits=9)    # Ensure that we have a fix length of digits (9 is for Kepler, TESS or PLATO need more)
        ProductsStarID=[]
        for star in ProductsRot["starID"]:
            ProductsStarID.append(format_ID(extract_ID_from_processname(star), Ndigits=9))
    else:
        ProductsStarID=ProductsRot["starID"]
    # Scan the list of starID to find all the indexes with matches
    indexesID=[]
    for i in range(len(ProductsStarID)):
        if ProductsStarID[i] == ID:
            indexesID.append(i)
    if indexesID == []:
        err_msg="Error: Cannot find the requested ID={} in ProductsRot['starID']".format(ID)
        raise ValueError(colored(err_msg, "red"))
    # Reduce the selection by looking at the ModelCode. 
    IndexUnique=[]
    for ind in indexesID:
        if ProductsRot["model"][ind] == ModelCode:
            IndexUnique.append(ind)
    if len(IndexUnique) == 1: 
        data=ProductsRot["data"][confidence][IndexUnique[0],:]
    else:
        print(colored("List of available ModelCode for ID={}".format(ID), "red"))
        for i in range(len(ProductsRot["model"])):
            print(colored("\t" + ProductsRot["model"][i], "red"))
        if len(IndexUnique) >1:
            err_msg="Error: Too many matches for ModelCode={} and ID={} in ProductsRot['model']\n. Debug Required.".format(ModelCode, ID)
            raise ValueError(colored(err_msg, "red"))
        elif IndexUnique ==[]:
            err_msg="Error: Cannot find the requested ModelCode={} in ProductsRot['model']".format(ModelCode)
            raise ValueError(colored(err_msg, "red"))
        else:
            err_msg="UNKNOWN Error in query_params_from_ProductsRot(). Debug required"
            raise ValueError(colored(err_msg, "red"))
    # Deal with ParamName input...
    if ParamName != None:
        try:
            ind_label=ProductsRot["label"].index(ParamName)
            data=data[ind_label]
        except ValueError:
            err_msg="Error: Cannot find the request parameter P={} for ModelCode={}, ID={}".format(ParamName, ModelCode, ID)    
            print(colored("Error: List of available parameters in ProductsRot['label']:", "red"))
            for l in ProductsRot["label"]:
                print("\t" + colored(l, "red"))
            raise ValueError(colored(err_msg,"red"))
    return data

def query_proba_from_ProductsOdds(ProductsOdds, ID, ModelCode=None, Confidence=None, ResolveID=True):
    '''
        A function that recover a specific array of probability associated to a given ID number. Optionally, it can also retrieve only a
        single probability value given a ModelCode
        ProductsOdds: The dictionary containing all of the data regarding Probabilities
        ID: The model identifier, either in full (in the form kplrXXXX) or in part (eg. only the KIC number)
        ModelCode: The model code associated to the model. Can be either 1001, 1101, 1111, 1201Gatedecompose-1 etc...
        Confidence: Define the name of the keys for the confidence. Default is all of them (usually "1s", "median", "1p")
        ResolveID (Optional): If set to True (default) will keep only the ID number and remove any other parts from the ID and from within the table 
        Returns: If Confidence is not specified, returns a dictionary with all the Confidence as keys. Otherwise, returns 
                 whose dimension depends on the number of ModelCode (dim=1 if ModelCode provided as argument)
    '''
    if ResolveID == True:
        ID=extract_ID_from_processname(ID)  # Ensure that ID is only numbers  
        ID=format_ID(ID, Ndigits=9)    # Ensure that we have a fix length of digits (9 is for Kepler, TESS or PLATO need more)
        ProductsStarID=[]
        for star in ProductsOdds["starID"]:
            ProductsStarID.append(format_ID(extract_ID_from_processname(star), Ndigits=9))
    else:
        ProductsStarID=ProductsOdds["starID"]

    if Confidence == None:
        Confidence=["m1s", "median", "p1s"]
    else:
        key=Confidence
    try:
        posID=ProductsStarID.index(ID)
    except ValueError:
        err_msg="Error: Cannot find the requested ID={} in ProductsOdds".format(ID)
        raise(colored(err_msg, "red"))

    if ModelCode != None:
        try:
            index=ProductsOdds["label"].index(ModelCode)
        except ValueError:
            err_msg="Error: Cannot find the requested ModelCode={} in ProductsOdds['label']".format(ModelCode)
            print(colored("List of available ModelCode in ProductsOdds['label']:", "red"))
            for l in ProductsOdds["label"]:
                print("\t" + colored(l, "red"))
            raise ValueError(colored(err_msg, "red"))
        if isinstance(Confidence, list) == True:
            data={}
            for key in Confidence:
                data[key]=ProductsOdds["Probabilities"][key][posID,index]
            data["label"]=ProductsOdds["label"]
        else:
            data=ProductsOdds["Probabilities"][key][posID,index]
    else: #    Case where all ModelCode is requested
        if isinstance(Confidence, list) == True: # Case where all the Confidences are requested
            data={}
            for key in Confidence:
                data[key]=ProductsOdds["Probabilities"][key][posID,:]
            data["label"]=ProductsOdds["label"]
        else:
            data=ProductsOdds["Probabilities"][key][posID,:]
    return data


def check_consistency(ProductsOdds, ProductsRot):
    '''
        This function just verify that for any 'starID' in ProductsOdds, we have a counterpart 'starID' in ProdcutsRot
        ProductsOdds : Dictionary read from Odds Ratio Product files generated by do_main_summary.py and organised by get_productsdata()
        ProductsRot : Dictionary read from Rotation Product files generated by do_main_summary.py and organised by get_productsdata()
        Return status: If True, the two dictionaries have a one-one relationship. Otherwise, something is broken
    '''
    status=True
    return status

def probability_composer(ProductsOdds, ProductsRot, Groups, DoRenormalise=False, EnforceNorm=True):
    '''
        This function takes the table of Probabilities obtained from the Odds ratio with the various models
        and recompose them into new more condensed classes of probabilities using the user-specified grouping of class. 
        For example, if the user want to regroup models with activity to answer the question: 
            Is there any activity at any latitude, independently from the underlying model?
        She/He would need to regroup models of type 1201* and 1211*
        ProductsOdds : Dictionary read from Odds Ratio Product files generated by do_main_summary.py and organised by get_productsdata()
        ProductsRot : Dictionary read from Rotation Product files generated by do_main_summary.py and organised by get_productsdata()
        Groups: List of lists containing all of the elements that must be regrouped. For example:
            [["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2"],
             ["1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]]
        DoRenormalise (Optional): If True, will renormalise the probabilities for all of the listed groups to 1. Otherwise, will keep the sum of the probabilities.
        EnforceNorm (Optional): If True, will raise an error if the sum of the probabilities is greater than 100 percent. 
                Otherwise, will raise a warning. This check is meaningful only for the Median values. Uncertainties on the probablity may make the value
                greater than 100 percent for key='m1s','m2s','p1s' or 'p2s'.
        Returns:
            RecomposedProductsOdds : A new dictionary that contains the Odds ratio, according to the regrouping rules
            RecomposedProductsRot : A new dictionary that may be altered to have average rotation parameters of the regrouped
                models (if do_average_ProductsRot = True). Or may not be changed (if do_average_ProductsRot = False)
    '''
    # Consistency check for the provided Group lists: Are all elements to regoup existing in the dictionaries? 
    for g in Groups:
        for model in g:
            found=False
            for l in ProductsOdds["label"]:
                if l == model:
                    found=True
            if found == False:
                err_msg="Error in probability_composer(): The model class {} providing the variable Group does not exist in ProductsOdds['label']".format(model)
                print(colored("The list of availabe models in ProductsOdds is:" , "red"))
                for l in ProductsOdds["label"]:
                    print("\t" + colored(l, "red"))
                raise ValueError(colored(err_msg,"red"))
            for l in ProductsRot["label"]:
                if l == model:
                    found=True
            if found == False:
                err_msg="Error in probability_composer(): The model class {} providing the variable Group does not exist in ProductsRot['label']".format(model)
                print(colored("The list of availabe models in ProductsRot is:" , "red"))
                for l in ProductsRot["label"]:
                    print("\t" + colored(l, "red"))
                raise ValueError(colored(err_msg,"red"))
    # Perform a hard copy of dictionaries if Rule 1 is applied (preserve array size). Otherwise, need to define a new container
    RecomposedProductsOdds = copy.deepcopy(ProductsOdds)
    RecomposedProductsRot = copy.deepcopy(ProductsRot)#

    keys=ProductsOdds["Probabilities"].keys()
    for key in keys:
        Pnew=[] # Used to store the new probabilities
        for g in Groups: # We travel within the list of lists here to perform the grouping for each sublist
            for i in range(ProductsOdds["Nrows"]): # We have to perform the grouping for each Row (each star) within the array
                Pn=0.
                for model in g: # Performing the Probability sum. Could perhaps be improved with numpy. But clearer like this
                    index=ProductsOdds["label"].index(model)
                    Pn=Pn + float(ProductsOdds["Probabilities"][key][i,index])/100. # Cast to float
                    #Pnorm[i]=Pnorm[i] + ProductsOdds["Probabilities"][key][i,index]
                Pnew.append(Pn)
                # Preserve the size of the array, but replace Probabilities with the Group Probability
                for model in g:
                    index=ProductsOdds["label"].index(model)
                    if DoRenormalise == True:
                        RecomposedProductsOdds["Probabilities"][key][i,index]=Pn # We keep it in fraction to enforce float computation
                    else:
                        RecomposedProductsOdds["Probabilities"][key][i,index]=Pn*100. # Here we convert it back to percentage. This is the final result
        if DoRenormalise == True:
            # Re normalise the probabilities to 1 using the new probabilities Pn calculated for each group
            Pnorm=np.zeros(ProductsOdds["Nrows"], dtype=float) # Used to renormalise the probabilities to 1 
            #print(colored("Pnorm=","light_cyan"), Pnorm, flush=True)
            for i in range(ProductsOdds["Nrows"]):
                indexes=[]
                # Compute the normalisation factor
                Pnew=[] # Used to store the new probabilities and compute Pnorm for each star
                for g in Groups:
                    for model in g:
                        index=ProductsOdds["label"].index(model)
                        indexes.append(index)
                        # Add to the normalisation factor only if the probability value was not added before
                        if RecomposedProductsOdds["Probabilities"][key][i,index] not in Pnew:
                            Pnorm[i]=Pnorm[i] + RecomposedProductsOdds["Probabilities"][key][i,index]
                            Pnew.append(RecomposedProductsOdds["Probabilities"][key][i,index])
                            #print(colored(" Index :", "light_cyan"), index, colored(" Pr:", "light_cyan"), RecomposedProductsOdds["Probabilities"][key][i,index], flush=True)
                #print(colored("Pnorm[{}] : ".format(i), "yellow"), Pnorm[i], flush=True)
                if Pnorm[i] >= 1.:
                    if EnforceNorm == True:
                        err_msg="Error: The sum of the probabilities for the median is greater than 100 percent. Cannot renormalise. Debug required."
                        for g in Groups:
                            for model in g:
                                index=ProductsOdds["label"].index(model)
                                print(colored("RecomposedProductsOdds['Probabilities'][{}][{},{}]={}".format(key,i,index, RecomposedProductsOdds["Probabilities"][key][i,index]), "red"))
                        raise ValueError(colored(err_msg, "red"))
                    else:
                        print(colored("Warning: The sum of the probabilities for the median is greater than 100 percent.", "yellow"), flush=True)
                        print("\t EnforceNorm is set to False. The Warning will be ignored.", flush=True)
                else:
                    # If there is no strange behaviour, we renormalise the probabilities
                    for index in indexes:                
                        #print(colored("\t Before renormalisation (should be a fraction): ", "yellow"), RecomposedProductsOdds["Probabilities"][key][i,index], flush=True)
                        if RecomposedProductsOdds["Probabilities"][key][i,index] != 0:
                            RecomposedProductsOdds["Probabilities"][key][i,index]=100.*RecomposedProductsOdds["Probabilities"][key][i,index]/Pnorm[i]
                        if np.isnan(RecomposedProductsOdds["Probabilities"][key][i,index]) == True:
                            err_msg=colored("Error: NaN value detected in RecomposedProductsOdds['Probabilities'][{}][{},{}]".format(key,i,index), "red")
                            raise ValueError(err_msg)
                        #print(colored("\t After renormalisation: (should be a percentage)", "blue"), RecomposedProductsOdds["Probabilities"][key][i,index], flush=True)
                    # Any Probability associated to labels that were not in one of the groups provided by the user will be set to 0
                    for index in range(len(RecomposedProductsOdds["label"])):
                        if index not in indexes:
                            RecomposedProductsOdds["Probabilities"][key][:, index] = 0
    return RecomposedProductsOdds, RecomposedProductsRot

def ProbabilityPeakGroup(ModelFamily, ProductsOdds, ProductsRot, ExternalData,
                       Ylabel, Xlabel, xlabel_plot, ylabel_plot, xvar_is_external, ExternalXerrlabel=None, 
                       ProbaThresholds=[50, 70, 90],
                       ColorsThresholds=["LightGray", "Gray", "Blue", "Green"],
                       DoRenormalise=[False,"1001"]):
    """
    Function that extract a model code 1201* or 1211* and then 
    creates a formated output using [] in order to show the highest probability in an encoded way.

    Args:
        ModelFamily (string): The Familiy of the model code. That is either 1201 or 1211.
        ProductsOdds (dict): The odds ratio products data. Defaults to None. If not provided, it will be loaded.
        ProductsRot (dict): The rotation products data. Defaults to None. If not provided, it will be loaded.
        xvar (str): The name of the x-variable.
        err_xvar (str): The name of the error in the x-variable.
        xvar_is_external (bool): Flag indicating if the x-variable is external.
        ProbaThresholds (list, optional): The probability thresholds for the odds ratio. Defaults to [50, 70, 90].
        ColorsThresholds (list, optional): The colors for the thresholds. Defaults to ["LightGray", "Gray", "Blue", "Green"].
    """
    # Imutable variables
    excluded_keys=['xlabel', 'ylabel', 'hline']
    ConfidenceKeys=["16", "50", "84"]
    # Define how we plot each subgroups
    Markers=["o", "s", "D", "o", "s", "D", "v", "^", "P","v", "^", "P"]
    Fills=['none', 'none', 'none', 'full', 'full', 'full', 'none', 'none', 'none', 'full', 'full', 'full']
    if ModelFamily not in ["1201", "1211", "12X1"]:
        raise ValueError("ModelFamily must be either 1201 or 1211.")
    elif ModelFamily == "1201":
        Groups=["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2",
                "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]
    elif ModelFamily == "1211":
        Groups=["1211Gatedecompose_-1", "1211Gatedecompose_1", "1211Gatedecompose_2",
                "1211Triangledecompose_-1", "1211Triangledecompose_1", "1211Triangledecompose_2"]
    else:
        Groups=["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2",
                "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2",
                "1211Gatedecompose_-1", "1211Gatedecompose_1", "1211Gatedecompose_2",
                "1211Triangledecompose_-1", "1211Triangledecompose_1", "1211Triangledecompose_2"]
    # Constructing all of the plot classes for the different ModelCodes within the provided Groups
    # This is constructed by NOT making subclasses by probability thresholds. Hence the ProbaThresholds=[100]
    iClass=0
    for ModelCode in Groups:
        if xvar_is_external == False:
            Class=DoClassesRotvsRot(ModelCode, ProductsOdds, ProductsRot, Ylabel, Xlabel,
                    xlabel_plot, ylabel_plot, ProbaThresholds=[100], 
                    ColorsThresholds=["Gray"], 
                    MarkerThresholds=["o"],FillThresholds=['none'],
                    do_hline=[True,0],  ConfidenceKeys=ConfidenceKeys)
        else:
            if ExternalXerrlabel == None:
                err_msg="Error: The external x-variable (xvar_is_external = True) is requested but no uncertainty column is provided. Debug required."
                raise ValueError(colored(err_msg, "red"))
            Class=DoClassesExternalvsRot(ModelCode, ProductsOdds, ProductsRot, ExternalData, Ylabel, Xlabel, ExternalXerrlabel,
                            xlabel_plot, ylabel_plot, ProbaThresholds=[100], 
                            ColorsThresholds=["Gray"], 
                            MarkerThresholds=["o"],FillThresholds=['none'],
                            do_hline=[True,0],  ConfidenceKeys=ConfidenceKeys)
        # Check that we indeed have only one key per class
        allowed_keys=[]
        for key in Class.keys():
            if key not in excluded_keys:
                allowed_keys.append(key)
        if len(allowed_keys) != 1:
            err_msg="Error: The class dictionary should contain only one key (key='A') in addition to the three 'global' keys. Debug required."
            print(colored("Class.keys()= ", "red"), Class.keys())
            raise ValueError(colored(err_msg, "red"))
        #ClassGroups.append(Class)
        # Identify the argmax of "Probability"
        # First, we initialise the Index_Startlist Array that contain all of the index for the highest probability
        # Then for each star, get the index of the position in ClassGroups that has the highest probability
        if iClass == 0:
            #Index_Starlist=np.zeros(len(Class["A"]["Probability"]), dtype=int)
            Pr_Starlist=np.zeros(len(Class["A"]["starID"]), dtype=float)
            Xall=np.zeros(len(Class["A"]["starID"]), dtype=float)
            Yall=np.zeros(len(Class["A"]["starID"]), dtype=float)
            Xerrall=np.zeros((2,len(Class["A"]["starID"])), dtype=float)
            Yerrall=np.zeros((2,len(Class["A"]["starID"])), dtype=float)
            # Each star must have a different Marker and Fill type because within each class, the ModelCode is not the same
            MarkerModelCode=[""]*len(Class["A"]["starID"])   #np.zeros(len(Class["A"]["Probability"]), dtype=str) 
            FillModelCode=[""]*len(Class["A"]["starID"])   #np.zeros(len(Class["A"]["Probability"]), dtype=str)
        # Update the final values if the Probability is higher only in the currently read "Class"
        for s in range(len(Class["A"]["starID"])):
            if Pr_Starlist[s] < Class["A"]["Probability"][s]:  
                Pr_Starlist[s]=Class["A"]["Probability"][s]
                if DoRenormalise[0] == True: # We renormalise the probabilities relative to the Reference probability list, if requested
                    Pr_ref=query_proba_from_ProductsOdds(ProductsOdds, Class["A"]["starID"][s], 
                            ModelCode=DoRenormalise[1], 
                            Confidence="median", 
                            ResolveID=True)
                    if Pr_ref > 0:
                        Odds=Pr_Starlist[s]/Pr_ref
                        Pr_Starlist[s]=100./(1. + 1./Odds)
                    else:
                        Pr_Starlist[s]=100.                    
                Xall[s]=Class["A"]["x"][s]
                Yall[s]=Class["A"]["y"][s]
                Xerrall[:,s]=Class["A"]["err_x"][:,s]
                Yerrall[:,s]=Class["A"]["err_y"][:,s]
                MarkerModelCode[s]=Markers[iClass]
                FillModelCode[s]=Fills[iClass]
        iClass=iClass+1
    # Now we have the Index_Starlist that contains the index of the highest probability for each star
    # We can now Initialise and fill the output class that keeps only the highest probability for each star
    do_hline=[True,0]
    classes={"xlabel": xlabel_plot,
            "ylabel": ylabel_plot,
            "hline": do_hline}
    # We make a list of unique ID from the ProductsRot["starID"] list
    starIDAll=[]
    for star in ProductsOdds["starID"]:
        s=format_ID(extract_ID_from_processname(star), Ndigits=9) # fix length using Kepler format
        if s not in starIDAll: # ensure uniqueness of the list
            starIDAll.append(s)
    key="A"
    for i in range(len(ProbaThresholds)+1): 
        if i==0:
            PosKeep=np.where(Pr_Starlist <= ProbaThresholds[i])[0]
            label="P <= {}".format(ProbaThresholds[i])
        elif i>0 and i<len(ProbaThresholds):
            PosKeep=np.where(np.bitwise_and(Pr_Starlist > ProbaThresholds[i-1], Pr_Starlist < ProbaThresholds[i]))[0]
            label="{} < P < {}".format(ProbaThresholds[i-1], ProbaThresholds[i])
        else:
            PosKeep=np.where(Pr_Starlist >= ProbaThresholds[i-1])[0]
            label="P >= {}".format(ProbaThresholds[i-1])
        classes[key]={"label": label,
                "color":ColorsThresholds[i],
                "marker":[MarkerModelCode[j] for j in PosKeep],
                "fillstyle": [FillModelCode[j] for j in PosKeep], 
                "positions": PosKeep,
                "x": Xall[PosKeep],
                "y": Yall[PosKeep],
                "err_x":Xerrall[:,PosKeep],
                "err_y":Yerrall[:,PosKeep],
                "Probability":Pr_Starlist[PosKeep],
                "starID":[starIDAll[j] for j in PosKeep]
            }
        key=increment_letter(key)
    return classes

def propag_err_theta_min_max(theta0, delta, err_theta0, err_delta):
    '''
    Use the error propagation law for gaussian errors in order to compute uncertainties 
    on theta_min and thetea_max coming from (theta0, delta)
    theta0: Value of the mean latitude of activity
    delta: Extension of the active band
    err_theta0: A 2-list with the minus and plus uncertainty on theta0
    err_delta : A 2-list with the minus and plus uncertainty on delta
    '''
    theta_max=theta0 + delta/2
    theta_min=theta0- delta/2
    err_theta_max= [np.sqrt(err_theta0[0]**2 + err_delta[0]**2/4),  np.sqrt(err_theta0[1]**2 + err_delta[1]**2/4)]
    err_theta_min= [np.sqrt(err_theta0[0]**2 + err_delta[0]**2/4),  np.sqrt(err_theta0[1]**2 + err_delta[1]**2/4)]
    return theta_min, theta_max, err_theta_min, err_theta_max

def data_Benomar2022():
    # epsilon, err_epsi_m, err_epsi_p, theta0, theta0_m, theta0_p, delta, delta_m, delta_p, theta_min, err_theta_min_m, err_theta_min_p, theta_max, err_theta_max_m, err_theta_max_p
    # WARNING: ERROR PROPAGATION WAS USED FOR THETA_MIN AND THETA_MAX COMPUTATION
    #          IT IS BETTER LATER TO USE THE FULL PDF TO ACCOUNT FOR CORRELATIONS
    # ACTIVE SUN:
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(76, 7, [7,8], [5,16])
    Active_Sun={'Teff':5777, 
                'a1':[421, 10, 10],
                'epsilon':[7.6*1e-4, 5.3*1e-4, 25.5*1e-4], 
                'theta0': [76, 7, 8], 
                'delta': [7, 5, 16], 
                'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
                'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]] 
                }
    # 16 Cyg A:
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(71, 4, [14,13], [3,15])
    #
    CygA={"Teff":5825, 
            'a1':[614, 37, 37],
            'epsilon':[5.3*1e-4, 4.1*1e-4, 19.6*1e-4], 
            'theta0': [71, 14, 13], 
            'delta':[4, 3, 15],  
            'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
            'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]]
            }
    # 16 Cyg B (Sol with largest probability and largest uncertainty):
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(67, 3, [25,14], [2,13])
    CygB_higha4={'Teff':5750, 
                    'a1':[607, 78, 78],
                    'epsilon':[4.6*1e-4, 3.7*1e-4, 20.8*1e-4], 
                    'theta0': [67, 25, 14], 
                    'delta':[3, 2, 13],  
                    'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
                    'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]] 
                    }
    # 16 Cyg B (Sol with lower probability and lowest uncertainty):
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(67, 3, [25,14], [2,13])
    CygB_lowa4={'Teff':0, 
                    'a1':[607, 78, 78],
                    'epsilon':[12.5*1e-4, 7.5*1e-4, 31.9*1e-4], 
                    'theta0': [58, 3, 3], 
                    'delta':[8, 6, 17],  
                    'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
                    'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]] 
                    }
    return Active_Sun, CygA, CygB_lowa4, CygB_higha4


def scoreclass2color(score, limits=[-1, 2.5, 3.5]):
    '''
        Receives the score-based table and in function of thresholds return a color:
        The classification is as follow:
            - if score < 0: Unclassified ==> RED
            - if 0 <= score <= 2.5: Bad ==> Light Gray
            - if 2.5 < score <= 3.5: Fair ==> Blue
            - if score > 3.5: Good  ==> Dark Green
        The score is expected to be on a 5-point basis
    '''
    if score <limits[0]:
        color='red'
    if 0 <= score <= limits[1]:
        color='lightgray'
    if limits[1] < score <= limits[2]:
        color='blue'
    if score > limits[2]:
        color='darkgreen'
    return color

def find_indices(lst, element):
    '''
    Returns all element of a list that match a given value or string
    lst: The list within which we search for values
    element: the value or string that is looked inside the list lst
    '''
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)


def main_a2CF_a2AR_group(dir_stars='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/products/a1a2a3a4_l0123_a1Gauss/', ids_file='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/id_list.txt', outfile='stars.summary.aj'):
    '''
    Function that measure a2CF and calculate a2AR=a2-a2CF for an ensemble of stars
    Then it writes the outputs into a structure file
    '''
    step_factor=3
    first_sample=0
    use_Anl_weight=True #False
    #
    found_dirs, found_id, missing_id=getdir_idlist(dir_stars, ids_file, ids_col=0)
    if len(missing_id) !=0:
        print('Warning: Some stars were not found in the root directory:')
        for m in missing_id:
            print(m)
        time.sleep(2)
    #
    i=0
    outputs='#{0:10} {1:10} {2:10} {3:10}'.format('id' ,  'a1'  ,  'err_a1_m',  'err_a1_p')
    outputs=outputs +'{0:10} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:10} {8:10}'.format('a2','err_a2_m','err_a2_p','a2_CF','err_a2_CF_m','err_a2_CF_p','a2_AR','err_AR_m','err_AR_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('a3','err_a3_m', 'err_a3_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('a4','err_a4_m', 'err_a4_p')
    outputs=outputs + '{0:10} {1:10} {2:10}\n'.format('inc','err_inc_m', 'err_inc_p')
    for d in found_dirs:
        a2_CF_all, a2_CF_l, a2_CF_mean12_samples, a2_AR_mean_samples=do_a2CF_a2AR(d, step_factor=step_factor, first_sample=first_sample, use_Anl_weight=use_Anl_weight)
        a2_CF_stats=make_stats(a2_CF_mean12_samples, confidence=[2.25,16,50,84,97.75])
        a2_AR_stats=make_stats(a2_AR_mean_samples, confidence=[2.25,16,50,84,97.75])
        a2_CF_err=make_error_from_stats(a2_CF_stats)
        a2_AR_err=make_error_from_stats(a2_AR_stats)
        aj_stats, inc_stats, aj_samples, inc_samples, Nsamples=get_aj_inc_star(d)
        aj_err=make_error_from_stats(aj_stats)
        inc_err=make_error_from_stats(inc_stats)
        #    id   a1    err_a1_m  err_a1_p
        outputs=outputs + '{0:10} {1:10.4f} {2:10.4f} {3:10.4f}'.format(found_id[i], aj_stats[0,2], aj_err[0,0], aj_err[1,0])
        #    a2  err_a2_m  err_a2_p     a2_CF   err_a2_CF_m   err_a2_CF_p    a2_AR  err_AR_m err_AR_p 
        outputs=outputs +  '{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f} {5:10.4f} {6:10.4f} {7:10.4f} {8:10.4f}'.format(aj_stats[1,2], aj_err[0,1], aj_err[1,1], a2_CF_stats[2], a2_CF_err[0], a2_CF_err[1], a2_AR_stats[2], a2_AR_err[0], a2_AR_err[1])
        #    a3  err_a3_m  err_a3_p 
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(aj_stats[2,2], aj_err[0,2], aj_err[1,2])
        #    a4  err_a4_m  err_a4_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(aj_stats[3,2], aj_err[0,3], aj_err[1,3])
        #    inc  err_inc_m  err_inc_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}\n'.format(inc_stats[2], inc_err[0], inc_err[1])
        i=i+1
    #
    f=open(outfile, 'w')
    f.write(outputs)
    f.close() 
    print("All data are saved in:", outfile)


def main_a2CF_a2AR():
    '''
    Function to quickly call  main_ajAlm_group in a proper manner
    '''
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/19992002_incfix_fast_Priorevalrange/'
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/20062009_incfix_fast/'
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/kplr012069424_kasoc-wpsd_slc_v1_a2a3a4_nol3/'
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/kplr012069449_kasoc-wpsd_slc_v1_a2a3a4_nol3/'
    dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/kplr012069449_kasoc-wpsd_slc_v1_a2a3a4_nol3_GU2/'
    step_factor=3
    first_sample=0
    use_Anl_weight=True #False
    a2_CF_all, a2_CF_l, a2_CF_mean12, a2_AR_mean=do_a2CF_a2AR(dir_core, step_factor=step_factor, first_sample=first_sample, use_Anl_weight=use_Anl_weight)
    #
    exit()

def plot_ajAlm(classes, file_out=None, ax=None):
    '''
        A function that focus only on ploting in a standardised way the aj or Alm as a function of M, R, Teff,etc...
        example: if x = M and y = a1, you get M(a1)
        ax: ploting zone
        x: x-axis
        y: y-axis
        file_out : File where the plot is written. Used only if ax is provided
        classes: A set of classes which contain how the information is going to be plotted. 
                 Typically, if A is the best class and C is the worse, you may end up with:
                    {   "hline": [True, 0],
                        "xlabel": "Teff (K)",
                        "ylabel": "a2 (nHz)",
                        "A": {
                            "label": "Good (P>90%)",
                            "starID": [],
                            "color":"Green",
                            "fillstyle": "full",
                            "marker": 'o', 
                            "positions:" [0,1,4,7],
                            "x":[],
                            "y":[],
                            "yerr": np.zeros((2,4)), # asymetrical errors
                            "xerr": np.zeros((2,4))  # symetrical errors
                            },
                        "B": {},
                        "C": {}
                    }
                With "label" : The name given within the legend
                     "starID" : The list of stars that match the crtiera for classes A,B,...,Z
                     "color" : The color of the symbols
                     "hline" : A list with first value telling if we plot an horizontal line and the second at which y-value we plot it
                     "xlabel": x-axis label
                     "ylabel": y-axis label
                     "positions": The positions associated with the "label". It must be of same size as the original x and y
                     "x": The filtered x axis
                     "y": The filtered y axis
                     "yerr" : A 2D array (2,N) with asymetrical errors. lower index is for lower bound. upper index is for upper bound
                     "xerr" : A 2D array (2,N) with asymetrical uncertainties
        ax: Ploting zone: If provided, the plot is not initialised and the function just defines how to plots elements 
    '''
    excluded_keys=['xlabel', 'ylabel', 'hline']
    allowed_keys=[]
    for key in classes.keys():
        if key not in excluded_keys:
            allowed_keys.append(key)
    if ax == None:
        fig_1d, ax = plt.subplots(1,1, figsize=(12, 6))

    for key in allowed_keys:
        # convert marker, color and fillstyle into arrays
        if isinstance(classes[key]["marker"], str):
            markers = [classes[key]["marker"]]*len(classes[key]["x"])
        else:
            markers=classes[key]["marker"]
        if isinstance(classes[key]["color"], str):
            colors = [classes[key]["color"]]*len(classes[key]["x"])
        else:
            colors=classes[key]["color"]
        if isinstance(classes[key]["fillstyle"], str):
            fillstyles = [classes[key]["fillstyle"]]*len(classes[key]["x"])
        else:
            fillstyles=classes[key]["fillstyle"]
        # plot each point individually with error bars
        if len(classes[key]["x"]) != 0:
            ax.errorbar(classes[key]["x"], classes[key]["y"], 
                    xerr=classes[key]["err_x"], yerr=classes[key]["err_y"], 
                    linestyle='', color=colors[0], marker="none", fillstyle="none",
                    label=classes[key]["label"])
            for i in range(len(classes[key]["x"])):
                ax.plot(classes[key]["x"][i], classes[key]["y"][i], 
                    marker=markers[i], 
                    fillstyle=fillstyles[i], 
                    color=colors[i], label=None)
            
        if (classes["hline"][0]):
            ax.axhline(classes["hline"][1], linestyle='--', color='black')
        ax.set_ylabel(classes["ylabel"])
        ax.set_xlabel(classes["xlabel"])
    if ax == None:
        # Handling legends
        ax.legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(file_out, dpi=300)
    #exit()

def increment_letter(letter, increment=1):
    '''
        A basic function that provide increment letters
        Can be used in a loop to generate alphabetical letter
    '''
    if letter == "Z" or letter == "z":
        raise ValueError("Error: Cannot increment beyond the letter z or Z.\n Provide an initial letter that is not the last letter of the alphabet")
    else:
        return chr(ord(letter) + increment)

def DoClassesRotvsRot(ModelCode, ProductsOdds, ProductsRot, Ylabel, Xlabel,
                      xlabel_plot, ylabel_plot, ProbaThresholds=[50, 75, 90], 
                      ColorsThresholds=["Gray", "dimgray", "Blue", "Green"], 
                      MarkerThresholds=["o", "o", "o", "o"],FillThresholds=['none','none', "full", "full"],
                      do_hline=[True,0],  ConfidenceKeys=["16", "50", "84"]):
    classes={"xlabel": xlabel_plot,
             "ylabel": ylabel_plot,
             "hline": do_hline}
    key="A"    
    # We make a list of unique ID from the ProductsRot["starID"] list
    starIDAll=[]
    for star in ProductsRot["starID"]:
        s=format_ID(extract_ID_from_processname(star), Ndigits=9) # fix length using Kepler format
        if s not in starIDAll: # ensure uniqueness of the list
            starIDAll.append(s)
    # Recover all of the stars in the table and the errors 
    Xall=np.zeros(len(starIDAll), dtype=float)
    Xerrall=np.zeros((2,len(starIDAll)), dtype=float)
    Yall=np.zeros(len(starIDAll), dtype=float)
    Yerrall=np.zeros((2, len(starIDAll)), dtype=float)
    Pall=np.zeros(len(starIDAll), dtype=float)    
    # Now we deal with the Probabilties and with Rotation parameters
    j=0
    for star in starIDAll:
        # Proba
        Pall[j]=query_proba_from_ProductsOdds(ProductsOdds, star, ModelCode=ModelCode, Confidence="median", ResolveID=True)
        # X-axis
        Xall[j]=query_params_from_ProductsRot(ProductsRot, star, ModelCode, ConfidenceKeys[1], ParamName=Xlabel, ResolveID=True)
        x=query_params_from_ProductsRot(ProductsRot, star, ModelCode, ConfidenceKeys[0], ParamName=Xlabel, ResolveID=True) 
        Xerrall[0,j]=Xall[j] - x # lower error bound
        x=query_params_from_ProductsRot(ProductsRot, star, ModelCode, ConfidenceKeys[2], ParamName=Xlabel, ResolveID=True) 
        Xerrall[1,j]=x - Xall[j] # upper error bound
        # Y-axis
        Yall[j]=query_params_from_ProductsRot(ProductsRot, star, ModelCode, ConfidenceKeys[1], ParamName=Ylabel, ResolveID=True)
        y=query_params_from_ProductsRot(ProductsRot, star, ModelCode, ConfidenceKeys[0], ParamName=Ylabel, ResolveID=True) 
        Yerrall[0,j]=Yall[j] - y # lower error bound
        y=query_params_from_ProductsRot(ProductsRot, star, ModelCode, ConfidenceKeys[2], ParamName=Ylabel, ResolveID=True) 
        Yerrall[1,j]=y - Yall[j] # upper error bound
        j=j+1
    # And Finally we deal with the color coding and the properties of the classes len(ProbaThresholds)+1
    # Classes corresponding to the ProbaThresholds
    if len(ProbaThresholds) > 1:
        for i in range(len(ProbaThresholds)+1): 
            if i==0:
                PosKeep=np.where(Pall <= ProbaThresholds[i])[0]
                label="P <= {}".format(ProbaThresholds[i])
            elif i>0 and i<len(ProbaThresholds):
                PosKeep=np.where(np.bitwise_and(Pall > ProbaThresholds[i-1], Pall < ProbaThresholds[i]))[0]
                label="{} < P < {}".format(ProbaThresholds[i-1], ProbaThresholds[i])
            else:
                PosKeep=np.where(Pall >= ProbaThresholds[i-1])[0]
                label="P >= {}".format(ProbaThresholds[i-1])
            classes[key]={"label": label,
                    "color":ColorsThresholds[i],
                    "marker": MarkerThresholds[i],
                    "fillstyle": FillThresholds[i],
                    "positions": PosKeep,
                    "x": Xall[PosKeep],
                    "y": Yall[PosKeep],
                    "err_x":Xerrall[:,PosKeep],
                    "err_y":Yerrall[:,PosKeep],
                    "Probability":Pall[PosKeep],
                    "starID":[starIDAll[j] for j in PosKeep]
                }
            key=increment_letter(key)
    else:
        print(colored("\t Warning: The number of Probability thresholds is only one. A single class will be created "), "yellow")
        PosKeep=np.where(Pall <= ProbaThresholds[0])[0]
        label="P <= {}".format(ProbaThresholds[0])
        classes[key]={"label": label,
                "color":ColorsThresholds[0],
                "marker": MarkerThresholds[0],
                "fillstyle": FillThresholds[0],
                "positions": PosKeep,
                "x": Xall[PosKeep],
                "y": Yall[PosKeep],
                "err_x":Xerrall[:,PosKeep],
                "err_y":Yerrall[:,PosKeep],
                "Probability":Pall[PosKeep],
                "starID":[starIDAll[j] for j in PosKeep]
            } 
    return classes


def DoClassesExternalvsRot(ModelCode, ProductsOdds, ProductsRot,ExternalData, RotYlabel, ExternalXlabel, ExternalXerrlabel,
                            xlabel_plot, ylabel_plot, ProbaThresholds=[50, 75, 90], 
                            ColorsThresholds=["Gray", "dimgray", "Blue", "Green"], 
                            MarkerThresholds=["o", "o", "o", "o"],FillThresholds=['none','none', "full", "full"],
                            do_hline=[True,0],  ConfidenceKeys=["16", "50", "84"]):
    '''
        A function that generate the classes according to the Probability threshold we provide and using ProductsRot and ExternalData.
        This is then used by show_ajAlm to generate classes and make appropriate color-coded plots
        ModelCode: This is the value of the model that are shown. For example if you compare "1001" and "1101", 
                   you may want to know if 1101 is more likely than 1001. Then you should set ModelName="1101" 
        ProductsOdds: The dictionary of all Odds ratio and Probabilities. Typically in the format: 
                    {"starID": {}, "Probabilities": {}, "header": None, "label": None, "Ncols":0., "Nrows":0.}
        ProductsRot: The dictionary of all Rotation information. Typically in the format: 
                    {"starID":{}, "data": {}, "header": None, "label": None, "unit": None, "Ncols":0., "Nrows":0} 
        ExternalData: A dictionary that contains the ExternalData["starID"], ExternalData["label"] and ExternalData["data"] for all Teff, M, R, etc...
        RotXlabel: A string specifying the label to be looked for within ProductsRot. This will make the y-axis
        ExternalYlabel: A string specifying the label to be looked for within ExternalData["label"]. This will make the x-axis
        ExternalYerrlabel: A string specifying the label for the getting the error column in ExternalData["label"]
        xlabel_plot: A string the will be defining the x-axis label on the plot
        ylabel_plot: A string the will be defining the y-axis label on the plot
        ProbaThresholds: A list of Probability thresholds that define the color for the ploted quantities
        ColorsThresholds: A list of Colors that will be used to plot the Probabilities. Must of be of size(ProbaThresholds) + 1
        do_hline: A list that specify if an horizontal line is plotted do_hline[0] and what y-value do_hline[1]
        ConfidenceKeys: A list of 3 values containing the confidence levels that are going to be shown. Default is 1sigma
    '''
    invalid_list=[""] # list of invalid values that will be checked within the ExternalData["data"] matrix of values
    classes={"xlabel": xlabel_plot,
             "ylabel": ylabel_plot,
             "hline": do_hline}
    key="A"
    # Looking up for the request parameter in the external table to make the x-axis
    Xl0=ExternalData["label"].index(ExternalXlabel) # Find the location of the specified Ylabel
    Xerrl0=ExternalData["label"].index(ExternalXerrlabel) # Find the location of the specified Yerrlabel
    # We make a list of unique ID from the ProductsRot["starID"] list
    starIDAll=[]
    for star in ProductsRot["starID"]:
        s=format_ID(extract_ID_from_processname(star), Ndigits=9) # fix length using Kepler format
        if s not in starIDAll: # ensure uniqueness of the list
            starIDAll.append(s)
    # Formating in a fix length the ID
    for i in range(len(ExternalData["starID"])):
        ExternalData["starID"][i]=format_ID(ExternalData["starID"][i])
    # Removing stars that are not in all of the tables or without values in the External table
    StarsSkipped=[]
    StarIDOK=[]
    IndexStarIDOK=[]
    IndexExternalStarIDOK=[]
    it=1
    for star in starIDAll:
        print("[{}/{}]  Looking for {}...".format(it, len(starIDAll), star))
        failed = False
        try:
            i0=ExternalData["starID"].index(star)
        except ValueError:
            failed=True
            StarsSkipped.append(star)
            print(colored("    >>> Warning: The star {} is in ProductsRot but not in ExternalData".format(star), "yellow"))
        if failed == False:
            value=ExternalData["data"][i0][Xl0]
            if value not in invalid_list:
                IndexExternalStarIDOK.append(i0)
                IndexStarIDOK.append(starIDAll.index(star))
                StarIDOK.append(star)
            else:
                StarsSkipped.append(star)
                print(colored("     >>> Warning: The ExternalData value for the parameter {} of the star {} is not a valid. The star will be skipped".format(value, star)))
        it=it+1
    # Recover all of the stars in the table and the errors 
    Xall=np.zeros(len(IndexStarIDOK), dtype=float)
    Xerrall=np.zeros(len(IndexStarIDOK), dtype=float)
    Yall=np.zeros(len(IndexStarIDOK), dtype=float)
    Yerrall=np.zeros((2, len(IndexStarIDOK)), dtype=float)
    Pall=np.zeros(len(IndexStarIDOK), dtype=float)
    j=0
    # Making the X-axis using ExternalData
    for i in IndexExternalStarIDOK:
        Xall[j]=float(ExternalData["data"][i][Xl0])
        Xerrall[j]=float(ExternalData["data"][i][Xerrl0])
        j=j+1
    # Now we deal with the Probabilties and with Rotation parameters
    j=0
    for i in IndexStarIDOK:
        Pall[j]=query_proba_from_ProductsOdds(ProductsOdds, starIDAll[i], ModelCode=ModelCode, Confidence="median", ResolveID=True)
        Yall[j]=query_params_from_ProductsRot(ProductsRot, starIDAll[i], ModelCode, ConfidenceKeys[1], ParamName=RotYlabel, ResolveID=True)
        y=query_params_from_ProductsRot(ProductsRot, starIDAll[i], ModelCode, ConfidenceKeys[0], ParamName=RotYlabel, ResolveID=True) 
        Yerrall[0,j]=Yall[j] - y # lower error bound
        y=query_params_from_ProductsRot(ProductsRot, starIDAll[i], ModelCode, ConfidenceKeys[2], ParamName=RotYlabel, ResolveID=True) 
        Yerrall[1,j]=y - Yall[j] # upper error bound
        j=j+1
    # And Finally we deal with the color coding and the properties of the classes len(ProbaThresholds)+1
    # Classes corresponding to the ProbaThresholds
    if len(ProbaThresholds) > 1:
        for i in range(len(ProbaThresholds)+1): 
            if i==0:
                PosKeep=np.where(Pall <= ProbaThresholds[i])[0]
                label="P <= {}".format(ProbaThresholds[i])
            elif i>0 and i<len(ProbaThresholds):
                PosKeep=np.where(np.bitwise_and(Pall > ProbaThresholds[i-1], Pall < ProbaThresholds[i]))[0]
                label="{} < P < {}".format(ProbaThresholds[i-1], ProbaThresholds[i])
            else:
                PosKeep=np.where(Pall >= ProbaThresholds[i-1])[0]
                label="P >= {}".format(ProbaThresholds[i-1])
            #print(" >>>> PosKeep : ", PosKeep)
            # Define symetrical errors for the ExternalData parameters as the external data always contain symetrical errors
            err_x=np.zeros((2,len(PosKeep))) 
            err_x[0,:]=Xerrall[PosKeep]
            err_x[1,:]=Xerrall[PosKeep]
            classes[key]={"label": label,
                    "color":ColorsThresholds[i],
                    "marker": MarkerThresholds[i],
                    "fillstyle": FillThresholds[i],
                    "positions": PosKeep,
                    "x": Xall[PosKeep],
                    "y": Yall[PosKeep],
                    "err_x":err_x,
                    "err_y":Yerrall[:,PosKeep],
                    "Probability":Pall[PosKeep],
                    "starID":[StarIDOK[j] for j in PosKeep]
                }
            key=increment_letter(key)
    else:
        print("\t Warning: The number of Probability thresholds is only one. A single class will be created")
        PosKeep=np.where(Pall <= ProbaThresholds[0])[0]
        label="P <= {}".format(ProbaThresholds[0])
        err_x=np.zeros((2,len(PosKeep))) 
        err_x[0,:]=Xerrall[PosKeep]
        err_x[1,:]=Xerrall[PosKeep]
        classes[key]={"label": label,
                "color":ColorsThresholds[0],
                "marker": MarkerThresholds[0],
                "fillstyle": FillThresholds[0],
                "positions": PosKeep,
                "x": Xall[PosKeep],
                "y": Yall[PosKeep],
                "err_x":err_x,
                "err_y":Yerrall[:,PosKeep],
                "Probability":Pall[PosKeep],
                "starID":[StarIDOK[j] for j in PosKeep]
            }        
    return classes

def CheckInstanceDataClass(classes, key, subkey):
    """
    CheckInstanceDataClass function that check the validity of the given subkey for a specific key in the classes dictionary.
    If valid, it returns a list of elements based on the subkey.
    This function is specifically used in the WriteDataClass() function to handle the color, marker, and fillstyle subkeys.
    Parameters:
    - classes: A dictionary containing the data classes.
    - key: The key in the classes dictionary.
    - subkey: The subkey to be checked.

    Returns:
    - elements: The validated elements based on the subkey.

    Raises:
    - ValueError: If the subkey is not 'color', 'marker', or 'fillstyle', or if the value associated with the subkey is not a string, list, or numpy array.
    """
    if subkey != "color" and subkey != "marker" and subkey != "fillstyle":
        err_msg = "The {} must be either 'color', 'marker' or 'fillstyle'".format(subkey)
        raise ValueError(colored(err_msg), "red")
    
    if isinstance(classes[key][subkey], str):
        elements = np.repeat(classes[key][subkey], len(classes[key]["starID"]))
    elif isinstance(classes[key][subkey], list) or isinstance(classes[key][subkey], np.ndarray):
        elements = classes[key][subkey]
    else:
        err_msg = "The {} must be a string, a list or an array".format(subkey)
        raise ValueError(colored(err_msg), "red")
    return elements

def WriteDataClass(fileout, classes, header="#Tabular representation of the data class used in the plots\n"):
    """
    Write the data class to a file in a tabular format.

    Args:
        fileout (str): The output file path.
        classes (dict): A dictionary containing the data class.
        header (str, optional): The header text to be written at the beginning of the file. Defaults to "#Tabular representation of the data class used in the plots\n".

    Returns:
        None
    """
    header=header + "#x= " + classes["xlabel"] +"\n"
    header=header + "#y= " + classes["ylabel"] +"\n"
    excluded_keys=['xlabel', 'ylabel', 'hline']
    allowed_keys=[]
    for key in classes.keys():
        if key not in excluded_keys:
            allowed_keys.append(key)
    variables=""
    for key in allowed_keys:
        colors=CheckInstanceDataClass(classes, key, "color")
        markers=CheckInstanceDataClass(classes, key, "marker")
        fillstyles=CheckInstanceDataClass(classes, key, "fillstyle")
        labels="!{0:<5} {1:<15} {2:<15} {3:<15} {4:<15} {5:<15} {6:<15} {7:<15} {8:<10} {9:<10} {10:<5} {11:<8} {12:<15}\n".format("key", "StarID", "x", "y", "err_x_inf","err_x_sup", "err_y_inf", "err_y_sup", "Pr", "color", "marker", "fillstyle", "label")
        for i in range(len(classes[key]["starID"])):
            variables=variables + " {0:<5} {1:<15} {2:<15.5f} {3:<15.5f} {4:<15.5f} {5:<15.5f} {6:<15.5f} {7:<15.5f} {8:<10.2f} {9:<10} {10:<8} {11:<10} {12:<}\n".format(key, classes[key]["starID"][i], 
                                                classes[key]["x"][i], classes[key]["y"][i], 
                                                classes[key]["err_x"][0][i],classes[key]["err_x"][1][i],
                                                classes[key]["err_y"][0][i],classes[key]["err_y"][1][i],
                                                classes[key]["Probability"][i],
                                                colors[i],
                                                markers[i],
                                                fillstyles[i],
                                                classes[key]["label"])
        f=open(fileout, "w")
        f.write(header)
        f.write(labels)
        f.write(variables)
        f.close()

def varname2label(varname):
    '''
        This function translates a set of variable names into a human readable label using the following dictionary and based 
        on the latex format compatible with matplotlib.
        This is mostly to be used with the show_ajAlm() function.
        Inputs:
            varname: A string that is the name of the variable
        Outputs:
            label: A string that is the human readable label using the latex notation for matplotlib
    '''
    if varname == "a1":
        label=r"$a_1$ (nHz)"
    elif varname == "a2":
        label=r"$a_2$ (nHz)"
    elif varname == "a3":
        label=r"$a_3$ (nHz)"
    elif varname == "a4":
        label=r"$a_4$ (nHz)"
    elif varname == "epsilon":
        label=r"$\epsilon$"
    elif varname == "theta0":
        label=r"$\theta_0$"
    elif varname == "delta":   
        label=r"$\delta$"
    elif varname == "theta_min":
        label=r"$\theta_{min}$"
    elif varname == "theta_max":
        label=r"$\theta_{max}$"
    elif varname == "Teff"  or varname == "Teff_SDSS":
        label=r"T$_{eff}$ (K)"
    elif varname == "M":
        label=r"M (M$_\odot$)"
    elif varname == "R":
        label=r"R (R$_\odot$)"
    elif varname == "inc": 
        label=r"Inc (°)"
    else:  
        label=varname
        print(colored("Warning: The variable name {} is not recognised".format(varname), "yellow"))
        print(colored("         We will use the variable name as label", "yellow"))
    return label

def Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir, 
                       var_name1, var_name2, xvar, err_xvar, xvar_is_external,
                       core_odds_files="Proba_summary", core_rot_files="rotation_stats",
                       odds_keys=["m1s", "median", "p1s"],confidence_rot=[2.25, 16, 50, 84, 97.75], ProbaThresholds=[50, 70, 90],
                       ProductsOdds=None, ProductsRot=None,
                       do_recompose=[False, [[],[]]]):
    """
    Plot the ajAlm data after loading the correct data sets. Do not get confused with plot_ajAlm() which is a function that only plots the data.

    Args:
        set_dir (str): The directory of the data set.
        ModelCode (str): The model code.
        dir_out (str): The output directory.
        imgfile (str): The name of the output image file.
        NonSeismic (bool): Flag indicating if non-seismic data is included.
        products_root_dir (str): The root directory of the products data.
        var_name1 (str): The name of the first variable.
        var_name2 (str): The name of the second variable.
        xvar (str): The name of the x-variable.
        err_xvar (str): The name of the error in the x-variable.
        xvar_is_external (bool): Flag indicating if the x-variable is external.
        core_odds_files (str, optional): The core name for odds ratio files. Defaults to "Proba_summary".
        core_rot_files (str, optional): The core name for rotation files. Defaults to "rotation_stats".
        odds_keys (list, optional): The keys for odds ratio elements present in product dictionaries. Defaults to ["m1s", "median", "p1s"].
        confidence_rot (list, optional): The confidence levels for rotation. Defaults to [2.25, 16, 50, 84, 97.75].
        ProbaThresholds (list, optional): The probability thresholds for the odds ratio. Defaults to [50, 70, 90].
        ProductsOdds (dict, optional): The odds ratio products data. Defaults to None. If not provided, it will be loaded.
        ProductsRot (dict, optional): The rotation products data. Defaults to None. If not provided, it will be loaded.
        do_recompose (list, optional): A list indicating if recomposition is done using the probability_composer() function. 
                Defaults to [False, [[],[]], False]. The first element is a flag that if set to True instruct the program to perform
                the renormalisation. The second is the list of lists that contain the groups. 
                For example:
                [["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2"],
                ["1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]].
                The third element is the flag that if set to True instruct the program to perform the renormalisation of the probabilities
                to 1 between the elements of the group.
                Note on the output: 
                    The output of the probability_composer() function are two dictionaries that will be of same shape as the 
                    original ProductsOdds and ProductsRot dictionaries. The difference is that each Probability element for the user-given
                    groups will be summed up. In the example given above, you threfore will have:
                        a. Pr("1201_Gate_decompose_-1") = Pr("1201_Gate_decompose_1") = Pr("1201_Gate_decompose_2")
                        b. Pr("1201_Triangle_decompose_-1") = Pr("1201_Triangle_decompose_1") = Pr("1201_Triangle_decompose_2")
                    If the renormalisation is done, the sum of all of the probabilities for all groups will be set to 1.
    """
    if ProductsOdds == None or ProductsRot == None:
        ProductsOdds, ProductsRot=get_productsdata(os.path.join(products_root_dir , set_dir), core_odds_files=core_odds_files, core_rot_files=core_rot_files,
                    odds_keys=odds_keys, confidence_rot=confidence_rot)
    if do_recompose[0] == True:
        ProductsOdds, ProductsRot=probability_composer(ProductsOdds, ProductsRot, do_recompose[1], 
                                                       DoRenormalise=do_recompose[2], EnforceNorm=False) # Use of EnforceNorm=False to avoid raising error for 'm1s','m2s','p1s' and 'p2s' keys

    if xvar_is_external == True:
        classes_var_name1=DoClassesExternalvsRot(ModelCode, ProductsOdds, ProductsRot,
                            NonSeismic, var_name1, xvar, err_xvar, varname2label(xvar), varname2label(var_name1), do_hline=[True,0],
                            ProbaThresholds=ProbaThresholds, 
                            ColorsThresholds=["Gray", "dimgray", "Blue", "Green"], 
                            MarkerThresholds=["o", "o", "o", "o"],FillThresholds=['none','none', "full", "full"])
        WriteDataClass(os.path.join(dir_out, imgfile + ".{}.res".format(var_name1)), classes_var_name1, header="#Tabular representation of the data class used in the plots\n#Data from {}\n#Comparing: {}\n".format(dir_out, set_dir))
        classes_var_name2=DoClassesExternalvsRot(ModelCode, ProductsOdds, ProductsRot,
                            NonSeismic, var_name2, xvar, err_xvar, varname2label(xvar), varname2label(var_name2), do_hline=[True,0],
                            ProbaThresholds=ProbaThresholds, 
                            ColorsThresholds=["Gray", "dimgray", "Blue", "Green"], 
                            MarkerThresholds=["o", "o", "o", "o"],FillThresholds=['none','none', "full", "full"])
        WriteDataClass(os.path.join(dir_out, imgfile + ".{}.res".format(var_name2)), classes_var_name2, header="#Tabular representation of the data class used in the plots\n#Data from {}\n#Comparing: {}\n".format(dir_out, set_dir))
    else:
        classes_var_name1=DoClassesRotvsRot(ModelCode, ProductsOdds, ProductsRot, var_name1, xvar,
                        varname2label(xvar), varname2label(var_name1), ProbaThresholds=ProbaThresholds, 
                        ColorsThresholds=["Gray", "dimgray", "Blue", "Green"], 
                        MarkerThresholds=["o", "o", "o", "o"],FillThresholds=['none','none', "full", "full"],
                        do_hline=[True,0])
        WriteDataClass(os.path.join(dir_out, imgfile + ".{}.res".format(var_name1)), classes_var_name1, header="#Tabular representation of the data class used in the plots\n#Data from {}\n#Comparing: {}\n".format(dir_out, set_dir))
        classes_var_name2=DoClassesRotvsRot(ModelCode, ProductsOdds, ProductsRot, 'a4', 'a1',
                        varname2label(xvar), varname2label(var_name2), ProbaThresholds=ProbaThresholds, 
                        ColorsThresholds=["Gray", "dimgray", "Blue", "Green"], 
                        MarkerThresholds=["o", "o", "o", "o"],FillThresholds=['none','none', "full", "full"],
                        do_hline=[True,0])
        WriteDataClass(os.path.join(dir_out, imgfile + ".{}.res".format(var_name2)), classes_var_name2, header="#Tabular representation of the data class used in the plots\n#Data from {}\n#Comparing: {}\n".format(dir_out, set_dir))

    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    plot_ajAlm(classes_var_name1, ax=ax[0])
    plot_ajAlm(classes_var_name2, ax=ax[1])
    ax[0].set_title("Significance {}".format(set_dir))
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper right')	
    fig_1d.savefig(os.path.join(dir_out, imgfile), dpi=300)

def Doplot_ajAlm_MaxProba12011211(set_dir, ModelFamilyActivity, ModelCodeRef, dir_out, imgfile, NonSeismic, products_root_dir, 
                       var_name1, var_name2, xvar, err_xvar, xvar_is_external,
                       core_odds_files="Proba_summary", core_rot_files="rotation_stats",
                       odds_keys=["m1s", "median", "p1s"],confidence_rot=[2.25, 16, 50, 84, 97.75], ProbaThresholds=[50, 70, 90],
                       ProductsOdds=None, ProductsRot=None):
    """
    Plot the ajAlm data showing only Max proba among a set 1201 or 1211, RELATIVE to another set (eg. 1001).
    This is performed after loading the correct data sets.
    Do not get confused with plot_ajAlm() which is a function that only plots the data.

    Args:
        set_dir (str): The directory of the data set.
        ModelFamilyActivity (str): The model family code for the analysis for which the max(Probability) is shown (either 1201 or 1211).
        ModelCodeRef (str): The model code for the analysis used as reference (eg. 1001).
        dir_out (str): The output directory.
        imgfile (str): The name of the output image file.
        NonSeismic (bool): Flag indicating if non-seismic data is included.
        products_root_dir (str): The root directory of the products data.
        var_name1 (str): The name of the first variable.
        var_name2 (str): The name of the second variable.
        xvar (str): The name of the x-variable.
        err_xvar (str): The name of the error in the x-variable.
        xvar_is_external (bool): Flag indicating if the x-variable is external.
        core_odds_files (str, optional): The core name for odds ratio files. Defaults to "Proba_summary".
        core_rot_files (str, optional): The core name for rotation files. Defaults to "rotation_stats".
        odds_keys (list, optional): The keys for odds ratio elements present in product dictionaries. Defaults to ["m1s", "median", "p1s"].
        confidence_rot (list, optional): The confidence levels for rotation. Defaults to [2.25, 16, 50, 84, 97.75].
        ProbaThresholds (list, optional): The probability thresholds for the odds ratio. Defaults to [50, 70, 90].
        ProductsOdds (dict, optional): The odds ratio products data. Defaults to None. If not provided, it will be loaded.
        ProductsRot (dict, optional): The rotation products data. Defaults to None. If not provided, it will be loaded.
    """
    if ProductsOdds == None or ProductsRot == None:
        ProductsOdds, ProductsRot=get_productsdata(os.path.join(products_root_dir , set_dir), core_odds_files=core_odds_files, core_rot_files=core_rot_files,
                    odds_keys=odds_keys, confidence_rot=confidence_rot)

    classes_var_name1=ProbabilityPeakGroup(ModelFamilyActivity, ProductsOdds, ProductsRot, NonSeismic,
                        var_name1, xvar, varname2label(xvar), varname2label(var_name1), xvar_is_external, ExternalXerrlabel=err_xvar,
                        ProbaThresholds=ProbaThresholds,
                        ColorsThresholds=["LightGray", "Gray", "Blue", "Green"], DoRenormalise=[True, ModelCodeRef])
    
    classes_var_name2=ProbabilityPeakGroup(ModelFamilyActivity, ProductsOdds, ProductsRot, NonSeismic,
                        var_name2, xvar, varname2label(xvar), varname2label(var_name1), xvar_is_external, ExternalXerrlabel=err_xvar,
                        ProbaThresholds=ProbaThresholds,
                        ColorsThresholds=["LightGray", "Gray", "Blue", "Green"], DoRenormalise=[True, ModelCodeRef])
    # Saving the data in files
    WriteDataClass(os.path.join(dir_out, imgfile + ".PeakGroup.{}.res".format(var_name1)), classes_var_name1, header="#Tabular representation of the data class used in the plots\n#Data from {}\n#Comparing: {}\n".format(dir_out, set_dir))
    WriteDataClass(os.path.join(dir_out, imgfile + ".PeakGroup.{}.res".format(var_name2)), classes_var_name2, header="#Tabular representation of the data class used in the plots\n#Data from {}\n#Comparing: {}\n".format(dir_out, set_dir))
    #print(colored("Data saved in files.", "green"))
    # Actual plot
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    plot_ajAlm(classes_var_name1, ax=ax[0])
    plot_ajAlm(classes_var_name2, ax=ax[1])
    ax[0].set_title("Significance {}".format(set_dir))
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper right')	
    fig_1d.savefig(os.path.join(dir_out, imgfile), dpi=300)

def show_ajAlm(products_root_dir, external_file, core_odds_files="Proba_summary", core_rot_files="rotation_stats",
                odds_keys=["m1s", "median", "p1s"],
                confidence_rot=[2.25, 16, 50, 84, 97.75],
                dir_out=''):
    '''
    Function that plots the values of a2 and a4 in function of
        - a1
        - Teff
        - M
        - R
    If an activity_summary_file (with epsilon, theta0, delta) is also provided, shows its quantities the same quantities listed above
    '''
    products_root_dir="/Users/obenomar/Work/dev/ajAlm/data/HighLevelProducts/Tcoef_1.7_1.86.77/statistical_summary/"
    external_file="/Users/obenomar/Work/dev/ajAlm/External_data/composite_table_Legacy_Kamiaka.csv"
    dir_out="/Users/obenomar/Work/dev/ajAlm/data/HighLevelResults/Tcoef_1.7_1.86.77/"
 
    ProbaThresholds=[50, 70, 90]

    print('Reading the file with M, R, Teff, etc... ', external_file)   
    NonSeismic=get_externaldata(external_file) # Matricial format converted into a dictionary for some information

    '''
    REMAINING ISSUE HERE: DATA BASE OF 1001 AND 1101 ARE NOT THE SAME AS THE ONE IN ALL : Probabilities are SOMETIMES VERY DIFFERENT
    # ================= TEST ZONE =================
    set_dir = "1001_vs_1101"
    ModelCode="1101"
    imgfile='a2significance_vs_a2a4_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "Teff_SDSS", "Tot_eTeff", True)
    # --
    imgfile='a2significance_vs_a2a4_a1.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "a1", "", False)

    set_dir = "All"
    ModelCode="1101" # This is the Model that we will show for the rotation values
    imgfile='Activitysignificance_vs_1101_Teff.jpg'
    renormalize=True # If True, the probability is renormalized to 1 between the models listed in do_recompose
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1101"]], renormalize])
    imgfile='Activitysignificance_vs_1101_1111_Teff.jpg'
    renormalize=True # If True, the probability is renormalized to 1 between the models listed in do_recompose
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1101", "1111"]], renormalize])
    '''
    #Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
    #             do_recompose=[True, [["1001"],
    #                                  ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])

    renormalize=True # If True, the probability is renormalized to 1 between the models listed in do_recompose
    set_dir = "All"
    ModelCode="1201Gatedecompose_-1"
    imgfile='Activitysignificance_vs_1201Gate-1_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", 
                                       "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])
    ModelCode="1201Gatedecompose_1"
    imgfile='Activitysignificance_vs_1201Gate1_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", 
                                       "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])

    ModelCode="1201Gatedecompose_2"
    imgfile='Activitysignificance_vs_1201Gate2_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", 
                                       "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])

    ModelCode="1201Triangledecompose_-1"
    imgfile='Activitysignificance_vs_1201Triangle-1_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", 
                                       "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])

    ModelCode="1201Triangledecompose_1"
    imgfile='Activitysignificance_vs_1201Triangle1_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", 
                                       "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])

    ModelCode="1201Triangledecompose_2"
    imgfile='Activitysignificance_vs_1201Triangle2_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", 
                                       "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])

    ModelFamilyActivity="1201"
    ModelCodeRef="1001"
    imgfile='HIGHESTActivitysignificance_vs_1001_Teff.jpg'
    Doplot_ajAlm_MaxProba12011211(set_dir, ModelFamilyActivity, ModelCodeRef, dir_out + "/1201vs1001/", imgfile, NonSeismic, products_root_dir,  
                                  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True)

    ModelFamilyActivity="1211"
    ModelCodeRef="1001"
    imgfile='HIGHESTActivitysignificance_vs_1001_Teff.jpg'
    Doplot_ajAlm_MaxProba12011211(set_dir, ModelFamilyActivity, ModelCodeRef, dir_out + "/1211vs1001/", imgfile, NonSeismic, products_root_dir,  
                                  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True)

    ModelFamilyActivity="1211"
    ModelCodeRef="1011"
    imgfile='HIGHESTActivitysignificance_vs_1011_Teff.jpg'
    Doplot_ajAlm_MaxProba12011211(set_dir, ModelFamilyActivity, ModelCodeRef, dir_out + "/1211vs1011/", imgfile, NonSeismic, products_root_dir,  
                                  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True)

    ModelFamilyActivity="12X1"
    ModelCodeRef="1001"
    imgfile='HIGHESTActivitysignificance_vs_1001_Teff.jpg'
    Doplot_ajAlm_MaxProba12011211(set_dir, ModelFamilyActivity, ModelCodeRef, dir_out + "/12X1vs1001/", imgfile, NonSeismic, products_root_dir,  
                                  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True)

    # ================= TEST ZONE (END) =================
    
    '''
    # ================= a2  significance using 1001 vs 1101 =================
    set_dir = "1001_vs_1101"
    ModelCode="1101"
    imgfile='a2significance_vs_a2a4_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "Teff_SDSS", "Tot_eTeff", True)
    # --
    imgfile='a2significance_vs_a2a4_a1.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "a1", "", False)

    # ==============  Activity significance that is Model Constrained : 1001 vs Sum_{j,k} 12j1k ===============
    set_dir = "All"
    ModelCode="1211Triangledecompose_1" # This is the Model that we will show for the rotation values
    imgfile='Activitysignificance_vs_1211_Teff.jpg'
    renormalize=True # If True, the probability is renormalized to 1 between the models listed in do_recompose
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "theta0", "delta", "Teff_SDSS", "Tot_eTeff", True,
                 do_recompose=[True, [["1001"],
                                      ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]], renormalize])
     
    # ==============   a3  significance using 1101 vs 1111 ===============
    set_dir = "1101_vs_1111"
    ModelCode="1111"
    imgfile='a3significance_vs_a2a4_Teff.jpg'
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "Teff_SDSS", "Tot_eTeff", True)
    # -- 
    imgfile="a3significance_vs_a2a4_a1.jpg"
    Doplot_ajAlm(set_dir, ModelCode, dir_out, imgfile, NonSeismic, products_root_dir,  "a2", "a4", "a1", "", False)
    '''

    '''
    #
    # a2AR and a4 function of M
    fig_1d.savefig(dir_out+ 'summary_a2a4_Mass.jpg', dpi=300)
    #
    # a2AR and a4 function of R
    fig_1d.savefig(dir_out+ 'summary_a2a4_Radius.jpg', dpi=300)
    #
    # a2AR and a4 function of Teff
    fig_1d.savefig(dir_out+ 'summary_a2a4_Teff.jpg', dpi=300)
    #
    # a2AR and a4 function of past a1
    fig_1d.savefig(dir_out+ 'summary_a2a4_pasta1.jpg', dpi=300)
   #
    # a2AR and a4 function of past a3
    fig_1d.savefig(dir_out+ 'summary_a2a4_pasta3.jpg', dpi=300)

    # a2AR and a4 function of new a1
    fig_1d.savefig(dir_out+ 'summary_a2a4_newa1.jpg', dpi=300)
   #
    # a2AR and a4 function of new a3
    fig_1d.savefig(dir_out+ 'summary_a2a4_newa3.jpg', dpi=300)

    # a2AR and a4 function of inclination
    fig_1d.savefig(dir_out+ 'summary_a2a4_inc.jpg', dpi=300)

    ###############################################################
    ###############################################################
    ########################### ajAlm #############################
    ###############################################################
    ###############################################################
    if activity_summary_file != None:
        # a2AR and a4 function of theta0 FROM ajAlm
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        #ax[0].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
        ax[0].errorbar(data_ajAlm[pos_fair, i_theta0], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(data_ajAlm[pos_good, i_theta0], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('a2AR (nHz)')
        ax[0].set_xlim(0, 90)
        #ax[1].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
        ax[1].errorbar(data_ajAlm[pos_fair, i_theta0], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(data_ajAlm[pos_good, i_theta0], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].axhline(0, linestyle='--', color='black')
        ax[1].set_ylabel('a4 (nHz)')
        ax[1].set_xlabel('theta0')
        ax[1].set_xlim(0, 90)
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_a2a4_theta0.jpg', dpi=300)

        # a2AR and a4 function of delta FROM ajAlm
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        #ax[0].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
        ax[0].errorbar(data_ajAlm[pos_fair, i_delta], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(data_ajAlm[pos_good, i_delta], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('a2AR (nHz)')
        ax[0].set_xlim(0, 45)
        #ax[1].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
        ax[1].errorbar(data_ajAlm[pos_fair, i_delta], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(data_ajAlm[pos_good, i_delta], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].axhline(0, linestyle='--', color='black')
        ax[1].set_ylabel('a4 (nHz)')
        ax[1].set_xlabel('delta')
        ax[1].set_xlim(0, 45)
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_a2a4_delta.jpg', dpi=300)

        # theta0 and delta and epsilon function of M
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(M_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(M_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[0].set_yscale('log')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[1].errorbar(M_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(M_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(M_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(M_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Mass')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_Mass.jpg', dpi=300)

        # theta0 and delta and epsilon function of R
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(R_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(R_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(R_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(R_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(R_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(R_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Radius')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_Radius.jpg', dpi=300)

        # theta0 and delta and epsilon function of Teff
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_Teff.jpg', dpi=300)

        # theta0 and delta and epsilon function of past a1
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(45, linestyle='--', color='black')
        ax[2].set_xlabel('Past a1 (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_pasta1.jpg', dpi=300)

        # theta0 and delta and epsilon function of past a3
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(a3past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a3past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(a3past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(a3past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(90, linestyle='--', color='black')
        ax[2].errorbar(a3past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(a3past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(45, linestyle='--', color='black')
        ax[2].set_xlabel('Past a3 (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_pasta3.jpg', dpi=300)

        # ------
        # Added on 14 Dec 2022
        # ------
        # theta_min and theta_max function of Teff
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
        ax[1].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta_min], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta_min+1],data_ajAlm[pos_fair_Alm, i_theta_min+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta_min], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta_min+1],data_ajAlm[pos_good_Alm, i_theta_min+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel(r'$\theta_{min}$')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta_max], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta_max+1],data_ajAlm[pos_fair_Alm, i_theta_max+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta_max], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta_max+1],data_ajAlm[pos_good_Alm, i_theta_max+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel(r'$\theta_{max}$')
        ax[2].set_ylim(0, 90)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetamin_thetamax_Teff.jpg', dpi=300)
        # theta_min , theta_max, theta0 in a single plot ONLY FAIR AND GOOD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['Teff'], CygA['Teff'], CygB_higha4['Teff']]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])
#
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
#        for p in pos_bad_Alm:
#            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+2] ]
#            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # Add the Sun, 16 Cyg A and B from Benomar 2022
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['Teff'], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_Teff_GoodFair.jpg', dpi=300)
        #
        # theta_min , theta_max, theta0 in a single plot WITH gray BAD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['Teff'], CygA['Teff'], CygB_higha4['Teff']]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])
#
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
        for p in pos_bad_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['Teff'], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_Teff_GoodFairBad.jpg', dpi=300)

        # ---------
        # ---------

        # theta_min , theta_max, theta0 in a single plot ONLY FAIR AND GOOD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['a1'][0], CygA['a1'][0], CygB_higha4['a1'][0]]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])

        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
#        for p in pos_bad_Alm:
#            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+2] ]
#            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # Add the Sun, 16 Cyg A and B from Benomar 2022
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['a1'][0], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel(r'$a_1$ (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_a1_GoodFair.jpg', dpi=300)
        #
        # theta_min , theta_max, theta0 in a single plot WITH gray BAD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['a1'][0], CygA['a1'][0], CygB_higha4['a1'][0]]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
        for p in pos_bad_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['a1'][0], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel(r'$a_1$ (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_a1_GoodFairBad.jpg', dpi=300)
    else:
        print('WARNING: No activity summary file provided. Plots skipped!')

    print('Image files saved in: ', dir_out)
    '''

'''
# Main part of the code
dir_core='/Users/obenomar/Work/tmp/paper2_asphericity/'
IDs_file=dir_core+ '/RealData/id_list.txt'
aj_summary_file=dir_core + '/RealData/stars.summary.ajfit'
#activity_summary_file='/Users/obenomar/Work/trash/stars.summary.ajAlm'
activity_summary_file='/Users/obenomar/Work/trash/stars.summary.14Dec2022.ajAlm'
activity_summary_file='/Users/obenomar/Work/trash/stars.summary.15Dec2022-MORESAMPLES.ajAlm'
aj_classification_file=dir_core + '/RealData/stars.classification.ajfit'
show_aj(dir_core, IDs_file, aj_summary_file, aj_classification_file, activity_summary_file=activity_summary_file)

print('All tasks done')
'''