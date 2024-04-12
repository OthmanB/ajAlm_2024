'''
	The functions here have for goal to separate the contribution
	on a2 of the centrifugal force from the activity
    This is useful in order to properly propagate uncertainties of a1 and Dnu into the a2(AR) term
    And therefore, activity inference just need rely on a2(AR), 
    effectively using a2(AR)=a2 - a2(CF) as observable instead of a2
'''
import numpy as np
from glob import glob
from pythonlibs.read_outputs_tamcmc import bin2txt
from pythonlibs.process_outputs_tamcmc_library import get_Dnu_samples
from pythonlibs.process_outputs_tamcmc_library import get_nu_samples,get_height_width_samples
from pythonlibs.acoefs import eval_acoefs
#from  scipy.io import readsav
#from pythonlibs.process_outputs_tamcmc_library import reduce_samples
#from pythonlibs.process_outputs_tamcmc_library import read_parameters_length_IDL

def Qlm(l,m):
	Dnl=2./3
	Qlm=(l*(l+1) - 3*m**2)/((2*l - 1)*(2*l + 3))
	return Qlm*Dnl

def eta0_fct(Dnu=None, rho=None, verbose=False):
	if rho == None:
		if Dnu !=None:
			Dnu_sun=135.1
			#numax_sun=3150.
			R_sun=6.96342e5 #in km
			M_sun=1.98855e30 #in kg
			rho_sun=M_sun*1e3/(4*np.pi*(R_sun*1e5)**3/3) #in g.cm-3
			rho=(Dnu/Dnu_sun)**2 * rho_sun
			if verbose == True:
				print('Computing rho using Dnu...')
				print(' rho =', rho)
		else:
			print('Error in eta0(): You need to provide at least one of the following argument:')
			print('                 - Dnu')
			print('                 - rho')
			exit()
	G=6.667e-8 # cm3.g-1.s-2
	#eta0=3./(4.*np.pi*rho*G) # second^2 # WRONG BECAUSE WE ASSUMED OMEGA ~ a1. IT SHOULD BE OMEGA ~ 2.pi.a1
	eta0=3.*np.pi/(rho*G)
	return eta0

def nu_CF(nu_nl, Dnu, a1, l, m, a1_unit='nHz'):
   eta0=eta0_fct(Dnu=Dnu, rho=None)
   if a1_unit == 'nHz':
         return eta0*nu_nl * (a1*1e-9)**2 * Qlm(l,m)
   if a1_unit == 'microHz':
         return eta0*nu_nl * (a1*1e-6)**2 * Qlm(l,m)
   if a1_unit != 'nHz' and a1_unit != 'microHz':
         print('a1 must be provided either in nHz or in microHz')
         print('use the a1_unit argument of the nu_CF() function to set it properly')
         exit()

def a2_CF(nu_nl, Dnu, a1, l):
	nu_nlm=[]
	for m in range(-l, l+1):
		perturb=nu_CF(nu_nl, Dnu, a1, l, m)
		nu_nlm.append(nu_nl + perturb)
	acoefs=eval_acoefs(l, nu_nlm)
	#print(nu_nlm)
	return acoefs[1] # returns only a2

def getdir_idlist(dir_stars, ids_file, ids_col=0):
    '''
        A function that look for the ID number listed in a file (IDs_file) that 
        are representing stars within a root directory (dir_stars)
        Stars with IDs that are not identified are flagged in a separate list than the one which are found
        dir_stars: The directory where the subdirectories for all stars are supposed to be
        IDs_file: File with all the kic numbers. It may contain several columns. 
                  Header and comments are allowed and must be marked using a # as first caracter
        IDs_col:  Defines the column that must contain the kic numbers. Default is 0.
    '''
    # Read the file to get a list of kic numbers
    found_id=[]
    found_dir=[]
    missing_id=[]
    with open(ids_file, 'r') as fh:
        # Skip initial comments that starts with #
        for line in fh:
            if not is_comment(line):
                #id_list.append(line.split()[ids_col]) # If it is not a commented line, split it and take the value in the kics_col coloun
                id=line.split()[ids_col]
                dir=glob(dir_stars + '/*' + id + '*')
                if len(dir) == 0:
                    missing_id.append(id)
                else:
                    found_id.append(id)
                    found_dir.append(dir[0])
                    if len(dir) > 1:
                        print('Warning: Multiple directory that match the requested ID found for ID: ', id)
                        print('The program will exit so that you can decide what to do')
                        print('Please be sure that you have only one directory for each ID')
                        exit()
    return found_dir, found_id, missing_id

def is_comment(line, comment_marker='#'):
	'''
		function to check if a line
	    starts with some character.
	    Here # for comment
	'''
	return line.startswith(comment_marker)


def get_a2_CF_samples(nu_nl_samples, a1_samples, dnu_samples, weight=None, verbose=False):
    '''
        Function that takes all of the samples for each mode and for a1 and derive from it
        The average a2(CF) term
        nu_nl_samples: A tensor of data containing frequencies and of size [0,lmax],[0:Nmodes],[0:Nsamples]
        a1_samples: samples for the a1 coefficient 
    '''
    Nsamples=len(nu_nl_samples[0,0,:])
    Nmodes=len(nu_nl_samples[0,:,0])
    lmax=len(nu_nl_samples[:,0,0])-1
    try:
        w_dimension = np.ndim(weight)
        if w_dimension == 2:
            w = np.zeros((lmax+1, Nmodes, Nsamples))
            w[:,:Nmodes,:] = weight[:,:,np.newaxis]
        elif w_dimension == 3:
            w = weight
        else:
            raise Exception('The weight array must be 2D or 3D.')
    except:
        w = np.ones((lmax+1, Nmodes, Nsamples))
    if verbose == True:
        print("         a. Configuration:")
        print("             - Nsamples    : ", Nsamples)
        print("             - Nmodes      : ", Nmodes)
        print("             - lmax        : ", lmax)
        print("             - weights     : ", weight)
    #
        print("         b. Derive a2(CF) for each modes and averaging over n...")
    a2_CF_l=np.zeros((lmax+1,Nsamples))
    a2_CF_all=np.zeros((lmax+1, Nmodes,Nsamples))
    a2_CF_mean12=np.zeros(Nsamples)
    s=0
    for s in range(Nsamples):
        Norm_nl=np.sum(w[1:3, :, s]) # Norm of l=1,2 only
        for l in range(1,lmax+1):
            Norm_l=np.sum(w[l,:,s])
            for n in range(Nmodes):
                r=a2_CF(nu_nl_samples[l,n,s], dnu_samples[s], a1_samples[s], l)
                a2_CF_l[l,s]=a2_CF_l[l,s]+w[l,n,s]*r/Norm_l   # Averaging over Nmodes the a2_CF coefficient. Uses weights w[l,n] defined by default to 1
                a2_CF_all[l,n,s]=r
                if l<=2:
                    a2_CF_mean12[s]=a2_CF_mean12[s] + w[l,n,s]*r/Norm_nl
    return a2_CF_all, a2_CF_l, a2_CF_mean12


def compute_a2_AR(a2, a2_CF_mean):
    '''
        Function that takes the average a2_CF and remove it from the a2 in order to evaluate a2_AR
    '''
    a2_AR_mean=a2-a2_CF_mean   
    return a2_AR_mean 


#def do_a2CF_a2AR_IDL(dir_core, step_factor=1, first_sample=0, use_Anl_weight=False):
    '''
        Main function that compute the average a2_CF and a2_AR from a serie of samples of a1, a2, nu_nl
        nu_nl is used to derive Dnu
        a1 approximates the average internal rotation
        All of the samples must have been acquired using the TAMCMC program and must have been
        saved into sav IDL files (using eg. IDLPostMCMC code)
        dir_core: directory on which the products of the IDLPostMCMC analysis can be found
        step_factor: Control how many samples are skipped within the ensemble of samples inside sav files. Greatly enhance the computation time is set to >1
        first_sample: Control the initial samples that is used for the analysis. Can be used to remove a potential Burn-in phase
        use_Anl_weight: If True, uses the median of the Amplitudes of the modes in order to weight the <a2_CF>_nl computation
    '''
'''
    dir_sav=dir_core + '/Files/'
    print("A. Processing star with the following inputs:")
    print("  - dir_core      : ", dir_core)
    print("  - dir_files     : ", dir_sav)
    print("  - step_factor   : ", step_factor)
    print("  - first_sample  : ", first_sample)
    print("  - use_Anl_weight: ", use_Anl_weight)
    #
    print("B. Extracting relevant information...")
    #
    print("   1. plength...")
    plength=read_parameters_length_IDL(dir_sav, 'plength.txt')
    i0_aj=sum(plength[0:6]) # First Position after Nf_el list
    print("   2. Extract and reduce samples of a1...(ASSUMES NO SLOPE: a1_1=0)")
    a1_samples, Nsize=read_sav(dir_sav, i0_aj)
    a1_samples=reduce_samples(a1_samples, step_factor, first_sample)
    a1_samples=a1_samples*1e3 # conversion to nHz
    print("   3. Extract and reduce samples of a2... (ASSUMES NO SLOPE: a2_1=0)")
    a2_samples, Nsize=read_sav(dir_sav, i0_aj+2)
    a2_samples=reduce_samples(a2_samples, step_factor, first_sample)
    a2_samples=a2_samples*1e3
    print("   4. Extract and reduce samples of nu(n,l)...")
    nu_nl_samples=get_nu_samples(dir_sav, plength, Nsize, verbose=True, step_factor=step_factor, first_sample=first_sample)
    print("   5. Extract A(n,l) from synthese files...")
    r=readsav(dir_core+'/synthese.sav')
    Anl=r['STAT_SYNTHESE_AMPLITUDE'][:,:,3] # Keep only the median
    if use_Anl_weight == True:
        Nmodes=len(nu_nl_samples[0,:,0])
        try:
            lmax=3
            weight=np.zeros((lmax+1, Nmodes))
            for l in range(lmax+1):
                for n in range(Nmodes):
                    weight[l,n]=Anl[n,l]
        except:
            lmax=2
            weight=np.zeros((lmax+1, Nmodes))
            for l in range(lmax+1):
                for n in range(Nmodes):
                    weight[l,n]=Anl[n,l]
        print("weight:", weight)
    else:
        weight=None
    #
    Nsamples=len(a1_samples)
    lmax=len(nu_nl_samples[:,0,0])-1
    print("C. Compute...")
    print("     1. Dnu...")
    dnu_samples=get_Dnu_samples(nu_nl_samples[0,:,:])
    print("       - Dnu = {0:0.3f} +/- {1:0.3f}".format(np.median(dnu_samples), np.std(dnu_samples)) )
    print("     2. a2_CF...")
    a2_CF_all, a2_CF_l, a2_CF_mean12=get_a2_CF_samples(nu_nl_samples, a1_samples, dnu_samples, weight=weight) # NOTE: a2 IS RETURNED IN MICROHZ
    a2_CF_all=a2_CF_all*1e3
    a2_CF_l=a2_CF_l*1e3
    a2_CF_mean12=a2_CF_mean12*1e3
    print("               <a2>nl             :  ~ {} +/- {} (nHz)".format(np.median(a2_samples), np.std(a2_samples))) 
    print(" -- ")
    for l in range(1, lmax+1): 
        print("             a2_CF(l={})                 :  ~ {} +/- {}  (nHz)".format(l, np.median(a2_CF_l[l,:]), np.std(a2_CF_l[l,:])))
    print("             <a2_CF>nl (l=1,2 only)          :  ~ {} +/- {} (nHz)".format(np.median(a2_CF_all[1:3,:,:]), np.std(a2_CF_all[1:3,:,:]))) 
    print("          <a2_CF>nl (l=1,2 only) (a2_mean12) :  ~ {} +/- {} (nHz)".format(np.median(a2_CF_mean12), np.std(a2_CF_mean12))) 
    print("             <a2_CF>nl (all l)               :  ~ {} +/- {} (nHz)".format(np.median(a2_CF_all[1:,:,:]), np.std(a2_CF_all[1:,:,:]))) 
    print(" -- ")
    print("     3. a2_AR...")
    a2_AR_mean=compute_a2_AR(a2_samples, a2_CF_mean12)
    print("             <a2_AR>nl             :  ~ {} +/- {} (nHz)".format(np.median(a2_AR_mean), np.std(a2_AR_mean))) 
    print(" Note: From all of these values, it is recommended to use a2_mean12 as it involves the proper error propagation")
    return a2_CF_all, a2_CF_l, a2_CF_mean12, a2_AR_mean
'''  

def do_a2CF_a2AR(dir_core, process_name, a1_samples, a2_samples, use_Anl_weight=False,
					phase="A", chain=0, first_index=0, last_index=-1, period=1,
					cpp_path="cpp_prg/", outtmpdir="tmp/", verbose=True, samples_unit="muHz"):
    '''
        Main function that compute the average a2_CF and a2_AR from a serie of samples of a1, a2, nu_nl
        nu_nl is used to derive Dnu
        a1 approximates the average internal rotation
        All of the samples are directly read using bin2txt
        dir_core: directory on which the mcmc data are
        step_factor: Control how many samples are skipped within the ensemble of samples. Greatly enhance the computation time is set to >1
        use_Anl_weight: If True, uses the median of the Amplitudes of the modes in order to weight the <a2_CF>_nl computation
        process_name:  Specify the name of the process, as it was defined in the config_presets.cfg file of the TAMCMC analysis. Default is None.
        phase (str, optional): Phase of the MCMC process. Default is "A".
        chain (int, optional): Parallel chain from which the MCMC samples are taken. Default is 0.
        first_index (int, optional): Index of the first sample. Can be used to remove burn-in. Default is 0.
        last_index (int, optional): Last index of the chain. Default is -1.
        period (int, optional):  Periodicity at which we pick samples. Can be used to reduce the size of the samples. Default is 1.
        cpp_path (str, optional):  Directory where the compiled bin2txt file of the TAMCMC program is located. Default is "cpp_prg/".
        outtmpdir (str, optional):  Temporary directory used by bin2txt to extract binary files into txt files. Note that files within it are discarded. Default is "tmp/".
        verbose (bool, optional): If True, show detailled information on each step. Otherwise, return results (almost) silently
    '''
    if samples_unit == "muHz":
         a1_samples=a1_samples*1e3
         a2_samples=a2_samples*1e3
    if verbose == True:
        print("A. Processing star with the following inputs:")
        print("  - dir_core      : ", dir_core)
        print("  - process_name : ", process_name)
    #    
        print("B. Extracting all mode information...")
        print("   1. Use bin2txt to extract all of the parameters...")
    smcmc, labels, isfixed, plength=bin2txt(dir_core, process_name, phase=phase, chain=chain, 
	    first_index=first_index, last_index=last_index, period=period, single_param_index=-1,
		erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir=outtmpdir, 
		get_plength=True, verbose=verbose)

    if verbose == True:    
        print("   2. Extract and reduce samples of nu(n,l)...")
    nu_nl_samples=get_nu_samples(smcmc, plength, verbose=verbose, step_factor=1, first_sample=0)
    if verbose == True:
        print("   3. Extract H(n,l), W(n,l) and A(n,l) from full set of samples...")
    h_samples, w_samples, A_samples=get_height_width_samples(smcmc, nu_nl_samples, plength, verbose=verbose, do_amplitude=True)

    if use_Anl_weight == True:
        weight=A_samples
    else:
        weight=None
    #
    Nsamples=len(a1_samples)
    lmax=len(nu_nl_samples[:,0,0])-1
    if verbose == True:
        print("C. Compute...")
        print("     1. Dnu...")
    dnu_samples=get_Dnu_samples(nu_nl_samples[0,:,:])
    if verbose == True:
        print("       - Dnu = {0:0.3f} +/- {1:0.3f}".format(np.median(dnu_samples), np.std(dnu_samples)) )
        print("     2. a2_CF...")
    a2_CF_all, a2_CF_l, a2_CF_mean12=get_a2_CF_samples(nu_nl_samples, a1_samples, dnu_samples, weight=weight) # NOTE: a2 IS RETURNED IN MICROHZ
    a2_CF_all=a2_CF_all*1e3
    a2_CF_l=a2_CF_l*1e3
    a2_CF_mean12=a2_CF_mean12*1e3
    print("               <a2>nl             :  ~ {} +/- {} (nHz)".format(np.median(a2_samples), np.std(a2_samples))) 
    print(" -- ")
    if verbose == True:
        for l in range(1, lmax+1): 
            print("             a2_CF(l={})                 :  ~ {} +/- {}  (nHz)".format(l, np.median(a2_CF_l[l,:]), np.std(a2_CF_l[l,:])))
    #
    if verbose == True:
        print("             <a2_CF>nl (l=1,2 only)          :  ~ {} +/- {} (nHz)".format(np.median(a2_CF_all[1:3,:,:]), np.std(a2_CF_all[1:3,:,:]))) 
    print("          <a2_CF>nl (l=1,2 only) (a2_mean12) :  ~ {} +/- {} (nHz)".format(np.median(a2_CF_mean12), np.std(a2_CF_mean12))) 
    if verbose == True:
        print("             <a2_CF>nl (all l)               :  ~ {} +/- {} (nHz)".format(np.median(a2_CF_all[1:,:,:]), np.std(a2_CF_all[1:,:,:]))) 
        print(" -- ")
        print("     3. a2_AR...")
    a2_AR_mean=compute_a2_AR(a2_samples, a2_CF_mean12)
    print("             <a2_AR>nl             :  ~ {} +/- {} (nHz)".format(np.median(a2_AR_mean), np.std(a2_AR_mean))) 
    if verbose == True:
        print(" Note: From all of these values, it is recommended to use a2_mean12 as it involves the proper error propagation")
    return a2_CF_all, a2_CF_l, a2_CF_mean12, a2_AR_mean
 

def show_version():
    print(" Update on 24 Feb 2023: ")
    print("    - Discontinuing IDL functions")
    print("    - Full integration of a2CF and a2AR computation using directly bin2txt")
    print(" Update on 12 Oct 2023: ")
    print("    - Renaming the functions using the sav files by adding '_IDL' at the end.")
    print("    - Creating functions that use directly the binary mcmc files")