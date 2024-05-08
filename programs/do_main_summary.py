# This directory regroups the routines performing statistical summaries
# on MCMC results. 
# The goal is to regroup on a single file per star, A table that contains rotational summaries along with 
# the evidence associated to them
# Dependences are:
#   - analyse_evidence.py  : Provides functions to compute evidences for all of the models
#   - show_pdf_aj_obs.py : Provides functions to compute the pdf of aj and ajAlm related quantities for all of the models
import os 
import shutil
import numpy as np
from termcolor import colored
from show_pdf_aj_obs import do_all_analyse_ajAlm
from analyse_evidence import compute_Odds_ratios, compile_results_comparison

def compare_all_evidences(data_dir, dir_out, do_list, file_summary, append="false", recompute_evidence=[False]):
    '''
        Function in charge of computing all of the odds ratio between a list of provided models types
        and given a phase and a root data directory
    '''
    txtsize="{:26s}"
    try:
        os.makedirs(dir_out)
    except FileExistsError:
        pass
    #print("      >>> Computing all Odds ratio...")
    modeltypes, oddsratios, oddsratios_err, starID=compute_Odds_ratios(data_dir, phase="A",
                                                            allowed_modeltype=do_list[0], 
                                                            allowed_modelsubtype1=do_list[1],
                                                            allowed_modelsubtype2=do_list[2],
                                                            recompute_evidence=recompute_evidence)
    # Create an horizontal label axis
    xlabels=txtsize.format("modelID")
    for modeltype in modeltypes:
        s_m=""
        for m in modeltype:
            s_m = s_m + str(m)
        xlabels=xlabels + txtsize.format(s_m)
    xlabels=xlabels +"\n"
    # Save the oddsratios and their errors in files
    Nmodels=len(oddsratios[0,0,:]) # squarred matrix of size Nmodels
    for s in range(len(starID)):
        file_out_med=os.path.join(dir_out, starID[s] + ".odds")
        file_out_err=os.path.join(dir_out, starID[s] + ".errodds")
        header="# Odds ratio computed by do_main_summary.py::compare_all_evidences() version {}\n".format(version(verbose=False))
        header=header + "# StarID={:<50s}\n".format(starID[s])
        # Generate the matrix
        string_med=""
        string_err=""
        for i in range(Nmodels):
            # Concatenate elements defining the model type
            ylabels=""
            for m in modeltypes[i]:
                ylabels = ylabels + str(m)
            # Compose the median matrix of the odds
            string_med=string_med + "{:<26s}".format(ylabels)
            for j in range(Nmodels):
                string_med=string_med + "{:<26.6f}".format(oddsratios[s, i,j])
            string_med=string_med + "\n"
            # Compose the error matrix of the odds
            string_err=string_err + "{:<26s}".format(ylabels)
            for j in range(Nmodels):
                string_err=string_err + "{:<26.6f}".format(oddsratios_err[s, i,j])
            string_err=string_err + "\n"
        # Write the matrices
        with open(os.path.join(dir_out, file_out_med), "w") as f:
            f.write(header)
            f.write(xlabels)
            f.write(string_med)
        with open(os.path.join(dir_out, file_out_err), "w") as f:
            f.write(header)
            f.write(xlabels)
            f.write(string_err)
    # This is to make final summaries for Probabilities in form of tables and images
    compile_results_comparison(starID, modeltypes, oddsratios, oddsratios_err, 
                               dir_out, file_summary=file_summary, append=append, withcolors=False)
    return modeltypes, oddsratios, oddsratios_err, starID


def do_products_database(data_dir, dir_out, 
                        data_sources=["Legacy", "Kamiaka2018"],
                        done_modeltype=[1001,1011,1101, 1111, 1201, 1211],
                        done_modelsubtype1=["Gate", "Triangle"],
                        done_modelsubtype2=["decompose_-1", "decompose_1", "decompose_2"],
                        cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.86.78_x86/",
                        recompute_evidence=[False], phase="A", period=20, tmpdir="../tmp/"):
    '''
        The main function that compute databases. Its purpose is to:
            (1) Compute all odds ratio and evidences and regroup them within single files
            (2) generate all of the pdfs and summary tables for rotation and inclination related parameters
    '''
    
    do_list=[]
    do_list.append(done_modeltype)
    do_list.append(done_modelsubtype1)
    do_list.append(done_modelsubtype2)
    file_odds="Proba_summary"
    append=False
    for d in data_sources:
        print("Processing stars labeled as {}...".format(d))
        print("     - Summary of the evidences and odds ratio...")
        diro=os.path.join(dir_out,'odds_ratio')
        print(colored("     - Computation of evidence and odds ratios outputs...", "green"))
        modeltypes, oddsratios, obsratios_err, starID=compare_all_evidences(data_dir+d, diro, 
                                                                do_list, file_odds, append=append, recompute_evidence=recompute_evidence)
        print(colored("     - Summary of the pdfs...","green"), flush=True)
        do_all_analyse_ajAlm(data_dir, dir_out, d, prefix="kplr", 
                            allowed_modeltype=done_modeltype, 
                            allowed_modelsubtype1=done_modelsubtype1, 
                            allowed_modelsubtype2=done_modelsubtype2, 
                            show_inc=True, 
                            cpp_path=cpp_path, phase="A", period=period,
                            append=append, tmpdir=tmpdir)
        append=True # All iteration i>1 are in append mode

def do_tests():
    '''
        Function that perform a test of whole do_products_database() function on a small subset of data
    '''
    cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.86.78_x86/"
    data_dir="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Outputs/Tcoef_1/"
    dir_out="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Products/Tcoef_1/statistical_summary/test/"
    #
    recompute_evidence=[False, True, cpp_path]
    done_modeltype=[1001, 1111, 1201, 1211]
    #done_modeltype=[1211]
    #done_modeltype=[1001, 1111]
    shutil.rmtree(dir_out)
    os.mkdir(dir_out)
    do_products_database(data_dir, dir_out, 
                        data_sources=["test_set", "test_set2"],
                        done_modeltype=done_modeltype,
                        done_modelsubtype1=["Gate", "Triangle"],
                        done_modelsubtype2=["decompose_-1", "decompose_1", "decompose_2"],
                        cpp_path=cpp_path,
                        recompute_evidence=recompute_evidence)
    print("Finished. Data saved in:", flush=True)
    print("   ", dir_out, flush=True)
    
def do_products_main():
    '''
        This is the final outputs generator. It call do_products_database() for all of the computed models
    '''
    cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.86.78_x86/"
    data_dir="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Outputs/Tcoef_1.7_1.86.77/"
    dir_out="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Products/Tcoef_1.7_1.86.77/statistical_summary/All_Evidence_1chain_new/"
    tmpdir="../tmp/"
    #
    #recompute_evidence=[True, True, cpp_path]
    recompute_evidence=[False, True, cpp_path]
    done_modeltype=[1001,1011, 1101, 1111, 1201, 1211]
    do_products_database(data_dir, dir_out, 
                        data_sources=["Legacy", "Kamiaka2018"],
                        done_modeltype=done_modeltype,
                        done_modelsubtype1=["Gate", "Triangle"],
                        done_modelsubtype2=["decompose_-1", "decompose_1", "decompose_2"],
                        cpp_path=cpp_path,
                        recompute_evidence=recompute_evidence, phase="A", period=10, tmpdir=tmpdir)
    
def version(verbose=True):
    v="1.10"
    if verbose == True:
        print("StatSummaryCompute version {}".format(v), flush=True)
        print("   - 1.10 (8 May 2024): Added writing theta_min and theta_max in the summary tables")
        print("   - 1.01: First version of the code")
    return v


do_products_main()
#do_tests()


