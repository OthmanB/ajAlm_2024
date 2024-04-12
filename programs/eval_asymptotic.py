'''
    Program that use linear fits of the asymptotic to extract the asymtptotic relation 
    of the p modes from the MCMC samples
    Note that the full asymptotic function is derived from the test function developed for the C++ code at 
        the C++ program at https://github.com/OthmanB/rescale_freq_pmodes.git
    This C++ program can also rescale frequencies, but this function is not used here (decomposition only)
'''

import os
import numpy as np
import sys
sys.path.append('pythonlibs')
from read_outputs_tamcmc import bin2txt
from process_outputs_tamcmc_library import get_Dnu_epsilon, decompose_nu_nl, get_nu_samples

def do_pdf_asymptotic(dir_tamcmc_outputs, process_name, phase='A', chain=0,
                      first_index=0, last_index=-1, period=10, single_param_index=-1,
                      cpp_path='bin/', get_d01=False, get_d02=False):
    verbose=False # Turn on for debug only
    outdir='tmp/'
    smcmc, labels, isfixed, plength=bin2txt(dir_tamcmc_outputs, process_name, phase=phase, chain=chain, 
                                first_index=first_index, last_index=last_index, period=period, single_param_index=single_param_index,
                                erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outdir, get_plength=True)
    # point towards frequencies only
    step_factor=1 # If one wants to reduce further the samples... but this can be achieved earlier with period>1
    nu_samples=get_nu_samples(smcmc, plength, verbose=verbose, step_factor=step_factor)
    if get_d01 == False and get_d02 == False:
        # determine by linear fits the asymptotic for l=0
        dnu_samples, epsilon_samples, nbest_samples=get_Dnu_epsilon(nu_samples[0,:,:])
        return dnu_samples, epsilon_samples, nbest_samples, [], [], [],[] # Last terms are to indicate that d01,d02, O2_l1, O2_l2 are not retrieved
    else:
        # determine by linear fits the asymptotic for l=0 + get d0l by decomposition
        # do this for each sample
        Nsamples=len(nu_samples[0,0,:])
        dnu_samples=np.zeros(Nsamples)
        epsilon_samples=np.zeros(Nsamples)
        nbest_samples=np.zeros(Nsamples)
        if get_d01 == True:
            d01_samples=np.zeros(Nsamples)
            O2_l1_samples=np.zeros(Nsamples)
            nbest_samples=np.zeros((2,Nsamples))
        else:
            d01_samples=[]
            O2_l1_samples=[]
        if get_d02 == True:
            d02_samples=np.zeros(Nsamples)
            O2_l2_samples=np.zeros(Nsamples)
            nbest_samples=np.zeros((2,Nsamples))
        else:
            d02_samples=[]
            O2_l2_samples=[]
            nbest_samples=[]
        for s in range(Nsamples):
            fl0=nu_samples[0,:,s] 
            if get_d01 == True:
                n_best, Dnu, epsilon, d0l, O2_term=decompose_nu_nl(1, fl0, fl=nu_samples[1,:,s], Cfactor=0.25, verbose=verbose)
                dnu_samples[s]=Dnu
                epsilon_samples[s]=epsilon
                d01_samples[s]=d0l
                nbest_samples[0,s]=n_best[0]
                O2_l1_samples[s]=np.sum(O2_term**2/len(O2_term)) # sum of residuals 
            if get_d02 == True:
                n_best, Dnu, epsilon, d0l, O2_term=decompose_nu_nl(2, fl0, fl=nu_samples[2,:,s], Cfactor=0.25, verbose=verbose)
                dnu_samples[s]=Dnu
                epsilon_samples[s]=epsilon
                d02_samples[s]=d0l
                nbest_samples[1,s]=n_best[0]
                O2_l2_samples[s]=np.sum(O2_term**2/len(O2_term)) # sum of residuals
        return dnu_samples, epsilon_samples, nbest_samples, d01_samples, d02_samples, O2_l1_samples, O2_l2_samples

def test_pdf_asymptotic():
    dir_tamcmc_outputs="../data_tests/outputs_1001/"
    process_name="kplr008006161_kasoc-psd_slc_v1_1001"
    phase="A"
    period=100
    cpp_path="bin/"
    print(" ----  Test of the simple retrieval of Dnu, epsilon ---")
    dnu_samples, epsilon_samples, nbest_samples, d01_samples, d02_samples, O2_l1_samples, O2_l2_samples=do_pdf_asymptotic(dir_tamcmc_outputs, process_name, phase=phase, chain=0,
                      first_index=0, last_index=-1, period=period, single_param_index=-1,
                      cpp_path=cpp_path, get_d01=False, get_d02=False)
    print("dnu_samples ", np.median(dnu_samples), np.std(dnu_samples))
    print("epsilon_samples ", np.median(epsilon_samples), np.std(epsilon_samples))
    print("nbest :", np.median(nbest_samples), np.std(nbest_samples))
    logic_test=all(not a for a in [d01_samples,d02_samples,O2_l1_samples,O2_l2_samples])
    print("Is d01, d02, O2_l1, O2_l2 empty?", logic_test)
    #
    print(" ----  Test of the full retrieval of Dnu, epsilon, d01, d02 ---")
    dnu_samples, epsilon_samples, nbest_samples, d01_samples, d02_samples, O2_l1_samples, O2_l2_samples=do_pdf_asymptotic(dir_tamcmc_outputs, process_name, phase=phase, chain=0,
                      first_index=0, last_index=-1, period=period, single_param_index=-1,
                      cpp_path=cpp_path, get_d01=True, get_d02=True)
    print("dnu_samples ", np.median(dnu_samples), np.std(dnu_samples))
    print("epsilon_samples ", np.median(epsilon_samples), np.std(epsilon_samples))
    print("nbest :", np.median(nbest_samples), np.std(nbest_samples))
    print("d01 :", np.median(d01_samples), np.std(d01_samples))
    print("d02 :", np.median(d02_samples), np.std(d02_samples))
    print("Average O2_l1 (residuals) :", np.median(O2_l1_samples), np.std(O2_l1_samples))
    print("Average O2_l2 (residuals) :", np.median(O2_l2_samples), np.std(O2_l2_samples))
    
def do_all_asymptotic(dir_tamcmc_outputs, summary_output_file, phase='A', period=200, chain=0, first_index=0, last_index=-1,
                      cpp_path='bin/', dir_filter='kplr'):
    '''
        Function that scans a specific path for all of its MCMC processes and then perform the 
        asymptotic determination and save everything in a summary file
    '''
    print("Scanning provided directory keeping ony those starting by {}...".format(dir_filter))
    process_names = [subdir for subdir in os.listdir(dir_tamcmc_outputs) if os.path.isdir(os.path.join(dir_tamcmc_outputs, subdir)) and subdir.startswith(dir_filter)]
    Nstars=len(process_names)
    print(" ----  Test of the full retrieval of Dnu, epsilon, d01, d02 ---")
    i=1
    txtforfile=''
    for process_name in process_names:
        print("  [{}/{}] {}...".format(i,Nstars,process_name))
        dnu_samples, epsilon_samples, nbest_samples, d01_samples, d02_samples, O2_l1_samples, O2_l2_samples=do_pdf_asymptotic(dir_tamcmc_outputs, process_name, phase=phase, chain=chain,
                      first_index=first_index, last_index=last_index, period=period, single_param_index=-1,
                      cpp_path=cpp_path, get_d01=True, get_d02=True)
        print("   >>> nbest_samples: {}  , {}".format(np.median(nbest_samples[0,:]), np.std(nbest_samples[0,:])) )
        txtforfile=txtforfile + '{0}  {1:14.3f}  {2:14.3f}  {3:14.6f}  {4:14.6f}  {5:14.1f}  {6:14.1f}  {7:14.5f}  {8:14.5f}  {9:14.5f}  {10:14.5f}  {11:14.5f}  {12:14.5f}  {13:14.5f}  {14:14.5f}\n'.format(
                                    process_name, 
                                    np.median(dnu_samples), np.std(dnu_samples), 
                                    np.median(epsilon_samples), np.std(epsilon_samples),
                                    np.median(nbest_samples[0,:]), np.std(nbest_samples[0,:]), 
                                    np.median(d01_samples), np.std(d01_samples), 
                                    np.median(d02_samples), np.std(d02_samples), 
                                    np.mean(O2_l1_samples), np.std(O2_l1_samples),
                                    np.mean(O2_l2_samples), np.std(O2_l2_samples)
                                    )
        i=i+1

    header='#Summary of the large separation, epsilon [0,1], d01~2D0, d02~6D0 and the residuals O2\n'
    header=header + '#Results are obtained by fitting the data with MCMC and then processing then using eval_asymptotic::do_all_asymptotic (python code)\n'
    header=header + '#The goodness of the determination is goodness=mean(Sum(O2^2)/Nmodes) while the stddev(goodness) defines the sensitivity of the goodness towards the fit parameters (curvature of the likelihood). A large goodness means an important curvature in the modes \n'
    header=header + '#Median are used as best estimator and the standard deviation on the samples gives the uncertainty\n'
    header=header + '#Note on best_n0: This is the best estimate for the first relative radial order. If the error is non-zero, that may indicate edges issues while determining epsilon [0,1]. Manual inspection may be then required\n'
    header=header + '#Configuration summary for do_all_asymptotic():\n'
    header=header + '#Scanned directory dir_tamcmc_outputs: {}\n'.format(dir_tamcmc_outputs)
    header=header + '#Retained phase: {}\n'.format(phase)
    header=header + '#Retained chain: {}\n'.format(chain)
    header=header + '#Resampling periodicity period: {}\n'.format(period)
    header=header + '#first_index : {}\n'.format(first_index)
    header=header + '#last_index : {}\n'.format(last_index)
    header=header + '#ID                                       Dnu           err(Dnu)       epsilon         err(epsilon)        best_n0_l0      err(best_n0_l0)     d01           err(d01)         d02            err(d02)      goodness_norm_l1  stdev(goodness_l1)  goodness_norm_l2 stdev(goodness_l2)\n'
    f=open(summary_output_file, 'w')
    f.write(header)
    f.write(txtforfile)
    f.close()    

# Legacy stars
#dir_tamcmc_outputs='../../../Kepler_FullScale_Analysis/Sun-like/Outputs/Legacy/Tcoef_1/1001/'
#summary_output_file='../../../Kepler_FullScale_Analysis/Sun-like/Outputs/Legacy/Tcoef_1/Summary_Legacy_1001.txt'
#do_all_asymptotic(dir_tamcmc_outputs, summary_output_file, phase='A', period=100, chain=0, first_index=0, last_index=-1,
#                      cpp_path='bin/', dir_filter='kplr')

# Kamiaka stars
dir_tamcmc_outputs='../../../Kepler_FullScale_Analysis/Sun-like/Outputs/Kamiaka2018/Tcoef_1/1001/'
summary_output_file='../../../Kepler_FullScale_Analysis/Sun-like/Outputs/Kamiaka2018/Tcoef_1/Summary_Kamiaka_1001.txt'
do_all_asymptotic(dir_tamcmc_outputs, summary_output_file, phase='A', period=100, chain=0, first_index=0, last_index=-1,
                      cpp_path='bin/', dir_filter='kplr')