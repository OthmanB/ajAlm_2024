# process outputs
# This file contains a serie of functions to calculate various mode properties
# such as the the large separation, getting the aj parameters or the activity in the ajAlm models

import numpy as np
from  scipy.io import readsav
import os
import os.path
from pythonlibs.read_outputs_tamcmc import bin2txt

def get_Dnu_samples(nu_nl0_samples):
    '''
        nu_nl0_samples: Matrix of data with all samples for all l=0 modes
    '''
    Nsamples=len(nu_nl0_samples[0,:])
    dnu_samples=np.zeros(Nsamples)
    for s in range(Nsamples):
        dnu_samples[s]=get_Dnu(nu_nl0_samples[:,s])
    return dnu_samples

def get_Dnu_epsilon(nu_nl0_samples):
    '''
        nu_nl0_samples: Matrix of data with all samples for all l=0 modes
    '''
    Nsamples=len(nu_nl0_samples[0,:])
    dnu_samples=np.zeros(Nsamples)
    epsilon_samples=np.zeros(Nsamples)
    n0_samples=np.zeros(Nsamples)
    for s in range(Nsamples):
        dnu_samples[s], epsilon_samples[s], n0_samples[s]=get_Dnu(nu_nl0_samples[:,s], getepsilon=True)
    return dnu_samples, epsilon_samples, n0_samples

def get_Dnu(nu_nl0, getepsilon=False):
    x=np.linspace(0, len(nu_nl0)-1, len(nu_nl0))
    e=np.polyfit(x, nu_nl0, 1) # First order polynomial fit
    if getepsilon == True:
        r=e[1]/e[0]
        n0=np.fix(r)
        epsi=r-n0
        return e[0],epsi, n0
    else:
        return e[0]

def get_nu_samples(samples, plength, verbose=False, step_factor=1, first_sample=0):
    '''
        If you have the full set of samples for all parameters, you can use this to extract 
        from it only the frequency ones
        samples: a 2D array (Nsamples, Nparams)
		plength: a 1D list or array with the positions for the different parameters
        verbose: If true, shows the median of the extracted frequency 
        step_factor: allows you to skip "step_factor" samples to reduce the final setsize
        first_sample: allows you to keep samples starting at "first_sample
    '''
    Nsamples=len(samples[:,0])
    lmax=plength[1]
    Nfl=plength[2:6]
    i0_l=[sum(plength[0:2])]
    for l in range(1,lmax+1):
        i0_l.append(sum(plength[0:2]) + sum(Nfl[0:l]))

    if step_factor == 1 and first_sample == 0:
        Newsize=Nsamples
    else:
        Newsize=int((Nsamples -first_sample)/step_factor);
    nu_nl_samples=np.zeros((lmax+1, max(Nfl), Newsize))
    for l in range(lmax+1):
        for n in range(Nfl[l]):
            fnl=samples[:,i0_l[l] + n]
            if Newsize == Nsamples:
                nu_nl_samples[l, n, :]=fnl
            else:
                nu_nl_samples[l,n,:]=reduce_samples(fnl, step_factor, first_sample)
            if verbose == True:
                print("    (l,n) = ({},{})  : ~ {}".format(l, n, np.median(fnl)))
    return nu_nl_samples

def get_height_width_l0_samples(samples, plength, verbose=False, step_factor=1, first_sample=0):
    '''
        If you have the full set of samples for all parameters, you can use this to extract 
        from it only the height and width of l=0 ones, for most of models
		Warning: This may not work with all models. Only models with the following assumptions will be compatible:
				- Heights are at the begining of the parameter set and only for l=0
				- l>0 heights are based on visibility ratio Vl. This exclude any RGB model based on the asymptotics
        samples: a 2D array (Nsamples, Nparams)
		plength: a 1D list or array with the positions for the different parameters
        verbose: If true, shows the median of the extracted frequency 
        step_factor: allows you to skip "step_factor" samples to reduce the final setsize
        first_sample: allows you to keep samples starting at "first_sample
    '''
    Nsamples=len(samples[:,0])
    Nh=plength[0]
    i0_h=0
    i0_w=np.sum(plength[0:7])
    
    if step_factor == 1 and first_sample == 0:
        Newsize=Nsamples
    else:
        Newsize=int((Nsamples -first_sample)/step_factor);
    h_nl_samples=np.zeros((Nh, Newsize))
    w_nl_samples=np.zeros((Nh, Newsize))
    for n in range(Nh):
        hnl=samples[:,i0_h + n]
        wnl=samples[:,i0_w + n]
        if Newsize == Nsamples:
            h_nl_samples[n, :]=hnl
            w_nl_samples[n,:]=wnl
        else:
            h_nl_samples[n,:]=reduce_samples(hnl, step_factor, first_sample)
            w_nl_samples[n,:]=reduce_samples(wnl, step_factor, first_sample)
        if verbose == True:
            print("    (l,n) = ({},{})  : H ~ {}   ,   W ~ {}".format(0, n, np.median(hnl), np.median(wnl)))
    return h_nl_samples, w_nl_samples

def get_height_width_samples(samples, nu_samples,plength, verbose=False, step_factor=1, first_sample=0, do_amplitude=True):
	'''
		Use of get_height_width_l0_samples() combined with Visibilities and interpolation at nu_samples to 
		get the (n,l) height and widths. 
		Warning: Works only on Regular Main Sequence models. Specifically the model_MS_Global_aj_HarveyLike()
	'''
	# General settings and initialisation
	Nsamples=len(samples[:,0])
	lmax=plength[1]
	Nfl=plength[2:6]
	if step_factor == 1 and first_sample == 0:
		Newsize=Nsamples
	else:
		Newsize=int((Nsamples -first_sample)/step_factor);
	height_nl_samples=np.zeros((lmax+1, max(Nfl), Newsize))
	width_nl_samples=np.zeros((lmax+1, max(Nfl), Newsize))
	if do_amplitude == True:
		A_nl_samples=np.zeros((lmax+1, max(Nfl), Newsize))
	Vl_samples=np.zeros((lmax+1, Newsize))
	Vl_samples[0,:]=1 # Vl0 = 1 by definition
	# Extract l=0 height and width
	h_n0_samples, w_n0_samples=get_height_width_l0_samples(samples, plength,
											verbose=False, step_factor=step_factor, first_sample=first_sample)
	height_nl_samples[0,:,:]=h_n0_samples
	width_nl_samples[0,:,:]=w_n0_samples
	if do_amplitude == True:
		A_nl_samples[0,:,:]=h_n0_samples*w_n0_samples
	# Extract the visibilities
	i0=plength[0]-1
	for l in range(1, lmax+1):
		Vl=samples[:,i0 + l] 
		if Newsize == Nsamples:
			Vl_samples[l, :]=Vl
		else:
			Vl_samples[l,:]=reduce_samples(Vl, step_factor, first_sample)
		if verbose == True:
			print("    V({})  : ~ {}".format(l, np.median(Vl)))
		for n in range(Nfl[l]):
			for j in range(Newsize):
				height_nl_samples[l, n, j]=np.interp(nu_samples[l,n,j], nu_samples[0,:,j],h_n0_samples[:,j])*Vl_samples[l,j]
				width_nl_samples[l, n, j]=np.interp( nu_samples[l,n,j], nu_samples[0,:,j],w_n0_samples[:,j])
			if do_amplitude == True:
				A_nl_samples[l,n,:]=height_nl_samples[l,n,:]*width_nl_samples[l,n,:]
				if verbose == True:
					print("    H({},{})  : ~ {}     W({},{})  : ~ {}     A({},{})  : ~ {}".format(l, n, np.median(height_nl_samples[l, n,:]),
                                                                                   				  l, n, np.median(width_nl_samples[l, n,:]), 
																								  l, n, np.median(A_nl_samples[l, n,:]) ))
			else:
				if verbose == True:
					print("    H({},{})  : ~ {}     W({},{})  : ~ {}".format(l, n, np.median(height_nl_samples[l, n,:]), 
                                                                             l, n, np.median(width_nl_samples[l, n,:])))
	if do_amplitude == True:
		return height_nl_samples, width_nl_samples, A_nl_samples
	else:
		return height_nl_samples, width_nl_samples


def reduce_samples(param, step_factor, first_sample):
    Nsamples=len(param)
    step=max([1, np.floor(step_factor)])
    Newsize=int((Nsamples -first_sample)/step_factor);
    reduced_param=np.zeros(Newsize)
    lcpt=0
    cpt=first_sample
    while lcpt<Newsize:
        reduced_param[lcpt]=param[cpt];
        cpt=cpt+step_factor;
        lcpt=lcpt+1;
    return reduced_param

def decompose_nu_nl(l, fl0, fl=[], Cfactor=0.25, verbose=False):
    '''
        TAKEN from https://github.com/OthmanB/rescale_freq_pmodes.git (which is a more comprehenside C++ code)
        Function that takes the p modes frequencies and decompose them
        ie, identifies, Dnu, epsilon, n, d0l and O(2) (higher order) terms
        This allows to extract the different contribution to the modes 
        and use the following properties:
            - Dnu and epsilon are determined by the l=0 frequencies
            - d0l + O(2) and the radial order n are determined by solving analytically the linear fit
                to the frequencies at the fixed pre-determined value of the slope = Dnu
            - O(2) is determined by assuming it to be a 0-mean perturbation around d0l
    '''
    if fl == []: # In that case, we decompose for the l=0
        fl=fl0
        l=0
    #
    # Compute Dnu and epsilon
    n_ref=np.linspace(0, len(fl0)-1, len(fl0))
    f=np.polyfit(n_ref, fl0, 1)
    Dnu=f[0]
    epsilon, n0_0=np.modf(f[1]/f[0])
    # Least square on l=1 frequencies with a fix slope
    #    1. Isolate d0l + O(2)
    func=fl - (epsilon + l/2)*Dnu # This is Dnu_ref.n + d0l + O(2)
    #    2. Determine a first guess for n
    e0, n0_l=np.modf(func[0]/Dnu)
    n_l=np.linspace(n0_l, n0_l + len(func)-1, len(func))
    #    3. Make a list of adjacent potential n and find the one that ensure |Dnu_ref| >> |d0l| + O(2) and d0l < 0
    n_all=np.zeros((4,len(n_l)))
    n_all[0,:]=n_l-np.repeat(1, len(n_l))
    n_all[1,:]=n_l
    n_all[2,:]=n_l+np.repeat(1,len(n_l))
    n_all[3,:]=n_l+np.repeat(2,len(n_l))
    sol=np.zeros(4)
    for n_s in range(4):
        sol[n_s]=np.mean(func) - Dnu*np.mean(n_all[n_s,:])
    if verbose == True:
        print("func = ", func)
        print("sol = ", sol)
    # EndFor debug
    n_s_optimum=np.where(np.abs(sol) < Dnu*Cfactor)[0] # Identify solutions close to 0 as |d0l|<<|Dnu_ref|
    n_best=n_all[n_s_optimum,:].flatten() #This is the list of identified n following the conditions in the where
    d0l=sol[n_s_optimum] # This is the best solution for d0l + O(2)
    if len(d0l) > 1:
        print('Error: multiple solutions of d0l found')
        print('       Debug required')
        print('      d0l =', d0l)
        exit()
    else:
        d0l=d0l[0]
    #    4. Identify the residual variation by assuming that O(2) is a random term of 0-mean
    O2_term=func- n_best*Dnu - np.repeat(d0l, len(n_best))
    if verbose == True:
        print('--')
        print('Best list of n matching conditions |Dnu_ref| >> |d0l| + O(2) and d0l < 0:')
        print(n_best)
        print('---')
        print('Identified best solution for d0l:')
        print(d0l)
        print('---')
        print('Residual of the constrained fit:')
        print(O2_term)
    #print("mean(O2):", np.mean(O2_term))
    #plt.plot(fl1, O2_term)
    #plt.show()
    #exit()
    return n_best, Dnu, epsilon, d0l, O2_term

def read_rot_aj(d, i0_aj, ignore_slope=False, idlfiles=False, 
                    process_name=None,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=1,
					cpp_path="cpp_prg/", outtmpdir="tmp/", 
					jmax=6):
	"""
	Read the rotation parameters of the MS_Global_aj_Harvey_like model (available in TAMCMC version>1.70)
	The use specifies the directory of the sav file (d), the initial index at which rotation parameters are stored
	(i0_aj) and can decide to skip the reading of the parameters that handle the slope with frequency of the aj parameters
	(ignore_slope) if they wish so. Note that grid analysis of a-coefficients show that although there is a slope, this one
	does not improve the inference on the activity region if aj ~ aj_CF + aj_AR(theta, delta) (See Benomar+2022). This is why
	by default the slope is usually fixed to 0 in the TAMCMC code.
	VERSION COMPATIBLE ONLY WITH BIN2TXT for cpp_version>1.85.0
	
	Parameters:
		d (str): Directory in which to look for the data. If idlfiles = True, it must be the directory with sav files. Otherwise, it must be the directory with the bin files (raw MCMC files)
		i0_aj (int): First index where the rotation parameters are stored
		ignore_slope (bool, optional): If True, does not retrieve the slope information on the aj parameters, which is often fixed to 0. Default is False.
		idlfiles (bool, optional): If True, load files stored in the IDL format created by the PostTAMCMC (gradually depreciated method). Otherwise, read directly the bin files (raw MCMC files). Default is False.
		process_name (str, optional): Required setting when idlfiles = False. Specify the name of the process, as it was defined in the config_presets.cfg file of the TAMCMC analysis. Default is None.
		phase (str, optional): Required setting when idlfiles = False. Phase of the MCMC process. Default is "A".
		chain (int, optional): Required setting when idlfiles = False. Parallel chain from which the MCMC samples are taken. Default is 0.
		first_index (int, optional): Required setting when idlfiles = False. Index of the first sample. Can be used to remove burn-in. Default is 0.
		last_index (int, optional): Required setting when idlfiles = False. Last index of the chain. Default is -1.
		period (int, optional): Required setting when idlfiles = False. Periodicity at which we pick samples. Can be used to reduce the size of the samples. Default is 1.
		cpp_path (str, optional): Required setting when idlfiles = False. Directory where the compiled bin2txt file of the TAMCMC program is located. Default is "cpp_prg/".
		outtmpdir (str, optional): Required setting when idlfiles = False. Temporary directory used by bin2txt to extract binary files into txt files. Note that files within it are discarded. Default is "tmp/".

	Returns:
		 numpy.ndarray: Samples of all of the aj parameters
    """
	#jmax=6 # Maximum order for the a-coefficient
	Naj=2 # Number of paramters for each aj
	if idlfiles == True:
		a1_param, Nsize=read_sav_IDL(d, i0_aj)
	else:
		samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
			first_index=first_index, last_index=last_index, period=period, single_param_index=i0_aj,
			erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
		Nsize=len(samples)
	if ignore_slope == False:
		aj_param=np.zeros((jmax, Naj, Nsize))
		cpt=i0_aj
		for j in range(0,jmax):
			for i in range(Naj):
				if idlfiles == True:
					samples, N=read_sav_IDL(d, i0_aj+cpt)
				else:
					samples, labels, isfixed, plength=bin2txt(d, process_name, phase=phase, chain=chain, 
							first_index=first_index, last_index=last_index, period=period, single_param_index=i0_aj+cpt,
							erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
				aj_param[j, i,:]=samples
				cpt=cpt+1
	else:
		aj_param=np.zeros((jmax, Nsize))
		for j in range(0,jmax):
				if idlfiles == True:
					samples, N=read_sav_IDL(d, i0_aj+Naj*j)
				else:
					samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
							first_index=first_index, last_index=last_index, period=period, single_param_index=i0_aj+Naj*j,
							erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
				aj_param[j,:]=samples[:,0]
	return aj_param



def read_rot_ajAlm(d, i0_aj, ignore_slope=False, idlfiles=False, 
                    process_name=None,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=1,
					cpp_path="cpp_prg/", outtmpdir="tmp/"):
	"""
	Read the rotation parameters of the MS_Global_ajAlm_Harvey_like model (available in TAMCMC version>1.85.0)

	Parameters:
		d (str): Directory in which to look for the data. If idlfiles = True, it must be the directory with sav files. Otherwise, it must be the directory with the bin files (raw MCMC files)
		i0_aj (int): First index where the rotation parameters are stored
		ignore_slope (bool, optional): If True, does not retrieve the slope information on the aj parameters, which is often fixed to 0. Default is False.
		idlfiles (bool, optional): If True, load files stored in the IDL format created by the PostTAMCMC (gradually depreciated method). Otherwise, read directly the bin files (raw MCMC files). Default is False.
		process_name (str, optional): Required setting when idlfiles = False. Specify the name of the process, as it was defined in the config_presets.cfg file of the TAMCMC analysis. Default is None.
		phase (str, optional): Required setting when idlfiles = False. Phase of the MCMC process. Default is "A".
		chain (int, optional): Required setting when idlfiles = False. Parallel chain from which the MCMC samples are taken. Default is 0.
		first_index (int, optional): Required setting when idlfiles = False. Index of the first sample. Can be used to remove burn-in. Default is 0.
		last_index (int, optional): Required setting when idlfiles = False. Last index of the chain. Default is -1.
		period (int, optional): Required setting when idlfiles = False. Periodicity at which we pick samples. Can be used to reduce the size of the samples. Default is 1.
		cpp_path (str, optional): Required setting when idlfiles = False. Directory where the compiled bin2txt file of the TAMCMC program is located. Default is "cpp_prg/".
		outtmpdir (str, optional): Required setting when idlfiles = False. Temporary directory used by bin2txt to extract binary files into txt files. Note that files within it are discarded. Default is "tmp/".

	Returns:
		 numpy.ndarray: Samples of all of the aj parameters
    """
	jmax=3 # Maximum number of aj coefficient. as we read a1,a3,a5 ==> jmax=3
	Naj=2 # Number of paramters for each aj
	i0_epsilon=i0_aj + jmax*Naj # epsilon is made of two parameters
	i0_theta=i0_aj + (jmax+1)*Naj 
	i0_delta=i0_aj + (jmax+1)*Naj + 1
	# get a1,a3,a5
	aj_param=read_rot_aj(d, i0_aj, ignore_slope=ignore_slope, idlfiles=idlfiles, 
					process_name=process_name,
					phase=phase, chain=chain,
					first_index=first_index, last_index=last_index, period=period,
					cpp_path=cpp_path, outtmpdir=outtmpdir, 
					jmax=jmax)
	Nsize=len(aj_param[0,:])
	# get epsilon
	if idlfiles == True:
		if ignore_slope == False:
			epsilon_samples=np.zeros((Naj,Nsize))
			samples, N=read_sav_IDL(d, i0_epsilon)
			epsilon_samples[0,:]=samples
			samples, N=read_sav_IDL(d, i0_epsilon+1)
			epsilon_samples[1,:]=samples
		else:
			epsilon_samples=np.zeros(Nsize)
			epsilon_samples, N=read_sav_IDL(d, i0_epsilon)	
	else:
		if ignore_slope == False:
			epsilon_samples=np.zeros((Naj,Nsize))
			samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
			    first_index=first_index, last_index=last_index, period=period, single_param_index=i0_epsilon,
				erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
			epsilon_samples[0,:]=samples
			samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
			    first_index=first_index, last_index=last_index, period=period, single_param_index=i0_epsilon+1,
				erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
			epsilon_samples[1,:]=samples
		else:
			epsilon_samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
			    first_index=first_index, last_index=last_index, period=period, single_param_index=i0_epsilon,
				erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
	# get theta0
	if idlfiles == True:
		theta_samples, N=read_sav_IDL(d, i0_theta)
	else:
		theta_samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
			    first_index=first_index, last_index=last_index, period=period, single_param_index=i0_theta,
				erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)	
	# get delta
	if idlfiles == True:
		delta_samples, N=read_sav_IDL(d, i0_delta)
	else:
		delta_samples, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
			    first_index=first_index, last_index=last_index, period=period, single_param_index=i0_delta,
				erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)	
	
	if ignore_slope == False:
		ajAlm_param=np.zeros(((jmax+1)+2, Naj, Nsize))
		ajAlm_param[0:jmax,:,:]=aj_param
		ajAlm_param[jmax:jmax+Naj,:,:]=epsilon_samples
		ajAlm_param[jmax+Naj,:,:]=theta_samples
		ajAlm_param[jmax+Naj+1,:,:]=delta_samples
	else:
		ajAlm_param=np.zeros(((jmax+1)+2, Nsize))
		ajAlm_param[0:jmax,:]=aj_param
		ajAlm_param[jmax,:]=epsilon_samples[:,0]
		ajAlm_param[jmax+1,:]=theta_samples[:,0]
		ajAlm_param[jmax+2,:]=delta_samples[:,0]
	return ajAlm_param


def get_aj_inc_star(rootdir, confidence=[2.25,16,50,84,97.75], ignore_slope=True,
					process_name=None,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=1,
					cpp_path="cpp_prg/", outtmpdir="tmp/"):
	""" 
	Analyzes the results of a single star and extracts information on aj coefficients.
	This information is useful for bias analysis or ensemble analysis of aj coefficients.
	This program reads directly the binary outputs of the TAMCMC code

	Parameters:
	rootdir (str): The root directory where the binary data files are located.
	confidence (list, optional): Confidence intervals for statistics calculation. Default is [2.25,16,50,84,97.75].
	ignore_slope (bool, optional): If True, does not retrieve the slope information on the aj parameters, which is often fixed to 0. Default is False.		
	idlfiles (bool, optional): If True, load files stored in the IDL format created by the PostTAMCMC (gradually depreciated method). Otherwise, read directly the bin files (raw MCMC files). Default is False.
	process_name (str, optional): Required setting when idlfiles = False. Specify the name of the process, as it was defined in the config_presets.cfg file of the TAMCMC analysis. Default is None.
	phase (str, optional): Required setting when idlfiles = False. Phase of the MCMC process. Default is "A".
	chain (int, optional): Required setting when idlfiles = False. Parallel chain from which the MCMC samples are taken. Default is 0.
	first_index (int, optional): Required setting when idlfiles = False. Index of the first sample. Can be used to remove burn-in. Default is 0.
	last_index (int, optional): Required setting when idlfiles = False. Last index of the chain. Default is -1.
	period (int, optional): Required setting when idlfiles = False. Periodicity at which we pick samples. Can be used to reduce the size of the samples. Default is 1.
	cpp_path (str, optional): Required setting when idlfiles = False. Directory where the compiled bin2txt file of the TAMCMC program is located. Default is "cpp_prg/".
	outtmpdir (str, optional): Required setting when idlfiles = False. Temporary directory used by bin2txt to extract binary files into txt files. Note that files within it are discarded. Default is "tmp/".

	Returns:
	tuple: A tuple containing the following elements:
		- aj_stats (numpy.ndarray): Statistics of aj coefficients at different confidence intervals.
		- inc_stats (numpy.ndarray): Statistics of inclination parameters.
		- aj_samples (numpy.ndarray): Samples of all aj coefficients.
		- inc_samples (numpy.ndarray): Samples of inclination parameters.
		- num_samples (int): Number of samples.
	"""
	# Getting all of the indexes required for the process
	print("		1.  Preparing data..")
	inc, labels, isfixed, plength=bin2txt(rootdir, process_name, phase=phase, chain=chain, 
			first_index=first_index, last_index=last_index, period=period, single_param_index=0,
			erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=True) 
	i0_aj=sum(plength[0:6]) # First Position after Nf_el list
	i0_inc=sum(plength[0:-2]) # The last parameter is before the extra parameters that are at the end ([-1] position)
	#
	# Get the rotation parameters in form of samples
	print("		2. Gathering rotation parameters...")
	aj_samples=read_rot_aj(rootdir, i0_aj, ignore_slope=ignore_slope, idlfiles=False, 
                    process_name=process_name,
					phase=phase, chain=chain,
					first_index=first_index, last_index=last_index, period=period,
					cpp_path=cpp_path, outtmpdir=outtmpdir) 
	jmax=len(aj_samples[:,0])
	#
	# Get the inclination in form of samples
	print("		3. Gathering inclination parameters...")
	inc_samples=read_inc(rootdir, i0_inc, idlfiles=False, 
                    process_name=process_name,
					phase=phase, chain=chain,
					first_index=first_index, last_index=last_index, period=period,
					cpp_path=cpp_path, outtmpdir=outtmpdir)
	#
	# Extract basic statistics from the aj samples and frequency samples
	print("		4. Get stats using the samples of all of the parameters...")
	Nconfidence=len(confidence)
	aj_stats=np.zeros((jmax, Nconfidence))
	for j in range(0, jmax):
		aj_stats[j,:]=np.percentile(aj_samples[j,:], confidence) #make_stats(aj_samples[j,:], confidence=confidence)
	inc_stats=np.percentile(inc_samples, confidence) #make_stats(inc_samples)
	#
	return aj_stats, inc_stats, aj_samples, inc_samples, len(aj_samples[0,:])

def get_ajAlm_inc_star(rootdir, confidence=[2.25,16,50,84,97.75], ignore_slope=True,
					process_name=None,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=1,
					cpp_path="cpp_prg/", outtmpdir="tmp/"):
	""" 
	Analyzes the results of a single star and extracts information on aj coefficients, activity model (epsilon,theta,delta) and the inclination.
	This program reads directly the binary outputs of the TAMCMC code

	Parameters:
	rootdir (str): The root directory where the binary data files are located.
	confidence (list, optional): Confidence intervals for statistics calculation. Default is [2.25,16,50,84,97.75].
	ignore_slope (bool, optional): If True, does not retrieve the slope information on the aj parameters, which is often fixed to 0. Default is False.		
	idlfiles (bool, optional): If True, load files stored in the IDL format created by the PostTAMCMC (gradually depreciated method). Otherwise, read directly the bin files (raw MCMC files). Default is False.
	process_name (str, optional): Required setting when idlfiles = False. Specify the name of the process, as it was defined in the config_presets.cfg file of the TAMCMC analysis. Default is None.
	phase (str, optional): Required setting when idlfiles = False. Phase of the MCMC process. Default is "A".
	chain (int, optional): Required setting when idlfiles = False. Parallel chain from which the MCMC samples are taken. Default is 0.
	first_index (int, optional): Required setting when idlfiles = False. Index of the first sample. Can be used to remove burn-in. Default is 0.
	last_index (int, optional): Required setting when idlfiles = False. Last index of the chain. Default is -1.
	period (int, optional): Required setting when idlfiles = False. Periodicity at which we pick samples. Can be used to reduce the size of the samples. Default is 1.
	cpp_path (str, optional): Required setting when idlfiles = False. Directory where the compiled bin2txt file of the TAMCMC program is located. Default is "cpp_prg/".
	outtmpdir (str, optional): Required setting when idlfiles = False. Temporary directory used by bin2txt to extract binary files into txt files. Note that files within it are discarded. Default is "tmp/".

	Returns:
	tuple: A tuple containing the following elements:
		- aj_stats (numpy.ndarray): Statistics of aj coefficients at different confidence intervals.
		- inc_stats (numpy.ndarray): Statistics of inclination parameters.
		- aj_samples (numpy.ndarray): Samples of all aj coefficients.
		- inc_samples (numpy.ndarray): Samples of inclination parameters.
		- num_samples (int): Number of samples.
	"""
	# Getting all of the indexes required for the process
	print("		1.  Preparing data..")
	inc, labels, isfixed, plength=bin2txt(rootdir, process_name, phase=phase, chain=chain, 
			first_index=first_index, last_index=last_index, period=period, single_param_index=0,
			erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=True) 
	i0_aj=sum(plength[0:6]) # First Position after Nf_el list
	i0_inc=sum(plength[0:-2]) # The last parameter is before the extra parameters that are at the end ([-1] position)
	#
	# Get the rotation parameters in form of samples
	print("		2. Gathering rotation parameters...")
	ajAlm_samples=read_rot_ajAlm(rootdir, i0_aj, ignore_slope=ignore_slope, idlfiles=False, 
                    process_name=process_name,
					phase=phase, chain=chain,
					first_index=first_index, last_index=last_index, period=period,
					cpp_path=cpp_path, outtmpdir=outtmpdir) 
	jmax=len(ajAlm_samples[:,0])
	#
	# Get the inclination in form of samples
	print("		3. Gathering inclination parameters...")
	inc_samples=read_inc(rootdir, i0_inc, idlfiles=False, 
                    process_name=process_name,
					phase=phase, chain=chain,
					first_index=first_index, last_index=last_index, period=period,
					cpp_path=cpp_path, outtmpdir=outtmpdir)
	#
	# Extract basic statistics from the aj samples and frequency samples
	print("		4. Get stats using the samples of all of the parameters...")
	Nconfidence=len(confidence)
	aj_stats=np.zeros((jmax, Nconfidence))
	for j in range(0, jmax):
		aj_stats[j,:]=np.percentile(ajAlm_samples[j,:], confidence) #make_stats(ajAlm_samples[j,:], confidence=confidence)
	inc_stats=np.percentile(inc_samples, confidence) #make_stats(inc_samples)
	#
	return aj_stats, inc_stats, ajAlm_samples, inc_samples, len(ajAlm_samples[0,:])


def read_inc(d, i0_inc, idlfiles=False,
					process_name=None,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=1,
					cpp_path="cpp_prg/", outtmpdir="tmp/"):
	'''
            Inputs:
                d: directory in which to look for the data. If idlfiles = True, it must be the directory with sav files. 
                    Otherwise, it must be the directory with the bin files (raw MCMC files)
                i0_inc: Index where the inclination parameter is stored
                ignore_slope: If True, does not retrieve the slope information on the aj parameters, which is often fixed to 0
                idlfiles: If True, load files stored in the IDL format created by the PostTAMCMC (gradually depreciated method)
                            Otherwise, read directly the bin files (raw MCMC files)
                process_name: Required setting when idlfiles = False. Specify the name of the process, as it was defined in the config_presets.cfg file of the TAMCMC analysis
                phase: Required setting when idlfiles = False. Phase of the MCMC process.
                chain: Required setting when idlfiles = False. Parallel chain from which the MCMC samples are taken
                first_index: Required setting when idlfiles = False. Index of the first sample. Can be used to remove burn-in
                last_index: Required setting when idlfiles = False. Last index of the chain.
                period: Required setting when idlfiles = False. Periodicity at which we pick samples. Can be used to reduce the size of the samples
                cpp_path: Required setting when idlfiles = False. Directory where the compiled bin2txt file of the TAMCMC program is located
                outtmpdir: Required setting when idlfiles = False. Temporary directory used by bin2txt to extract binary files into txt files. Note that files within it are discarded. 
            Outputs:
                samples of all of the aj parameters	
	'''
	if idlfiles == True:
		inc, Nsize=read_sav_IDL(d, i0_inc)
	else:
		inc, labels, isfixed=bin2txt(d, process_name, phase=phase, chain=chain, 
				first_index=first_index, last_index=last_index, period=period, single_param_index=i0_inc,
				erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, get_plength=False)
	return inc

# --------------------------------------------------------------
#       FUNCTIONS ONLY DESIGNED FOR READING POST-MCMC OUTPUTS 
#                      (gradually depreciated)
# --------------------------------------------------------------
def get_aj_inc_star_IDL(rootdir, confidence=[2.25,16,50,84,97.75]):
	""" 
	Analyzes the results of a single star and extracts information on aj coefficients.
	This information is useful for bias analysis or ensemble analysis of aj coefficients.
	This program requires a pre-processing of the outputs using the postIDL MCMC code.

	Parameters:
	rootdir (str): The root directory where the data files are located.
	confidence (list, optional): Confidence intervals for statistics calculation. Default is [2.25,16,50,84,97.75].

	Returns:
	tuple: A tuple containing the following elements:
		- aj_stats (numpy.ndarray): Statistics of aj coefficients at different confidence intervals.
		- inc_stats (numpy.ndarray): Statistics of inclination parameters.
		- aj_samples (numpy.ndarray): Samples of all aj coefficients.
		- inc_samples (numpy.ndarray): Samples of inclination parameters.
		- num_samples (int): Number of samples.
	"""
	# Getting all of the indexes required for the process
	print("		1.  Preparing data..")
	param_file='plength.txt'
	plength=read_parameters_length_IDL(rootdir + '/Files/', param_file) # Read the parameters_length file and retrieves plength
	#Nf_el=plength[2:6] # Elements 2,3,4 and 5
	#i0_freq=sum(plength[0:2]) # Sum of elements 0 and 1 which are Nmax and lmax
	#print("plength: ", plength)
	i0_aj=sum(plength[0:6]) # First Position after Nf_el list
	i0_inc=sum(plength[0:-2]) # The last parameter is before the extra parameters that are at the end ([-1] position)
	#
	# Get the rotation parameters in form of samples
	print("		2. Gathering rotation parameters...")
	aj_samples=read_rot_aj(rootdir + 'Files/', i0_aj, ignore_slope=True, idlfiles=True) 
	jmax=len(aj_samples[:,0])
	#
	# Get the inclination in form of samples
	print("		3. Gathering inclination parameters...")
	inc_samples=read_inc(rootdir + 'Files/', i0_inc)
	#
	# Extract basic statistics from the aj samples and frequency samples
	print("		4. Get stats using the samples of all of the parameters...")
	Nconfidence=len(confidence)
	aj_stats=np.zeros((jmax, Nconfidence))
	for j in range(0, jmax):
		aj_stats[j,:]=np.percentile(aj_samples[j,:], confidence)
	inc_stats=np.percentile(inc_samples, confidence) #make_stats(inc_samples)
	#
	return aj_stats, inc_stats, aj_samples, inc_samples, len(aj_samples[0,:])


def read_sav_IDL(dir, ind, extra=None):
	'''
		Read an IDL sav file that is generated by the IDL PostMCMC program
	'''
	s, err=gen_data_filename_IDL(dir, ind, extra=extra)
	if err == True:
		print(s)
		exit()
	else:
		r=readsav(s)
	return r['param'], len(r['param'])

def read_parameters_length_IDL(dir, param_file):
	'''
		Read the text file that contains the parameters structure and that is generated by the IDL PostMCMC program
	'''
	plength=[]
	with open(dir + param_file, 'r') as f:
		# Skip initial comments that starts with #
		for line in f:
				plength.append(line) # Each line contains a single value

	plength=np.array(plength, dtype=int)
	return plength

def get_dirs_list(rootdir):
	'''
		A tiny function that scan a directory in order to find all of 
		its sub directories
	'''
	#print('os.scandir(rootdir)=', os.scandir(rootdir))
	dirs = [f.path for f in os.scandir(rootdir) if f.is_dir()]
	return dirs

def gen_data_filename_IDL(d, index, extra=None):
	'''
		Generate a string with a syntax consistent with filenames generated by the IDL PostMCMC program 
	'''
	if extra == None:
		extra=''
	err=False
	if index < 10:
		s=d+'00' + str(index) + extra + '.sav'
	if index >= 10 and index < 100:
		s=d+'0' + str(index) + extra + '.sav'
	if index >=100 and index < 1000:
		s=d + str(index) + extra + '.sav'
	if index >=1000:
		s='Cannot handle index numbers greater than 999'
		err=True
	return s, err