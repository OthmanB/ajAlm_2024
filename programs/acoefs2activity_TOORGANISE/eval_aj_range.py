import numpy as np
from scipy import interpolate
from fit_a2sig import Qlm
from fit_a2sig import eta0_fct
from acoefs import eval_acoefs
from read_outputs_tamcmc import bin2txt
from process_outputs_tamcmc_library import get_nu_samples, get_Dnu_samples, uncertainties
from io_tune import get_dirs_list
from termcolor import colored
import os

def read_matrix_tab(file, delete0=False):
	'''
		Reader of tables in a matrix format with a header marked by '#' as first character
	'''
	f=open(file)
	raw=f.read()
	f.close()
	raw=raw.split('\n')
	Nlines=len(raw)
	header=[]
	i=0
	for line in raw:
		if line != '':
			if line[0] == '#':
				header.append(line)
			else: 
				# If i=0, we need to initialise the matrix of data. Otherwise, write the Matrix of data
				s=line.split()
				if i==0:
					data=np.zeros((Nlines, len(s)), dtype=float)  
				data[i,:]=np.asarray(s, dtype=float)				
				i=i+1
	# Remove potential empty lines at the end
	if i-1 < Nlines and delete0 == True:
		data=data[0:i]
	return data, header

def read_id_list(file):
	'''
		Read a list of kic organized in a single column within a file.
		Assumes a single number per line
	'''
	f=open(file)
	raw=f.read()
	f.close()
	id_list=raw.split('\n')
	return id_list

def match_id_within_lists(myID, data_bellinger, data_pinsonnault, data_benomar, verbose=True):
	'''
		A small program that take the full table of parameters given by Bellinger+2019, Pinsonnault+2012 and from my list
		and try to match it with an user-provided ID.
		If it fails, it return an error
	'''
	i_ID=0 # Where the ID number is supposed to be in the tables

	#pos_bellinger=list(np.asarray(data_bellinger[:,i_ID], dtype=int)).index(int(myID))
	pos_bellinger=np.where(np.asarray(data_bellinger[:,i_ID], dtype=int) == int(myID))
	pos_pinsonnault=np.where(np.asarray(data_pinsonnault[:,i_ID], dtype=int) == int(myID))
	pos_benomar=np.where(np.asarray(data_benomar[:,i_ID],dtype=int) == int(myID))
	print('Searching ID:', myID)
	err=False
	if len(pos_bellinger[0]) != 0:
		if verbose == True:
			print('       ID found in Bellinger+2019 at index', pos_bellinger)
	else:
		if verbose == True:
			print(colored('       Error while searching in the Bellinger+2019 Table A1: could not find the requested ID', 'red'))
		err=True
	if len(pos_pinsonnault[0]) != 0:
		if verbose == True:
			print('       ID found in Pinsonnault+2012 at index', pos_pinsonnault)
	else:
		if verbose == True:
			print(colored('       Error while searching in the Pinsonnault+2012 Table 7: could not find the requested ID', 'red'))
		err=True
	if len(pos_benomar[0]) != 0:
		if verbose == True:
			print('       ID found in Benomar+2018 at index', pos_benomar)
	else:
		if verbose == True:
			print(colored('       Error while searching in the Benomar+2018 My table of a1 and a3: could not find the requested ID', 'red'))
		err=True
	if err == True:
		i_id=[]
		label=[]
		#raise ValueError(" Error while trying to match star catalogs")
		#exit()
	else:
		i_id=[pos_bellinger[0], pos_pinsonnault[0], pos_benomar[0]]
		label=['Bellinger+2019 Table A1', 'Pinsonnault+2012 Table 7', 'Benomar+2018']
	return i_id, label, err
	
def load_Alm_grids(dir_grid, gridfiles, lmax=2):
	Alm=[]
	if lmax>=1:
		Alm_grid_l1=np.load(dir_grid + gridfiles[0])
		Alm.append(Alm_grid_l1['Alm'])# Note: ftype = 'gauss' or 'gate' depends on the Alm grid file. No need to specify it
	if lmax >=2:
		Alm_grid_l2=np.load(dir_grid + gridfiles[1])
		Alm.append(Alm_grid_l2['Alm'])# Note: ftype = 'gauss' or 'gate' depends on the Alm grid file. No need to specify it
	if lmax >=3:
		Alm_grid_l3=np.load(dir_grid + gridfiles[2])
		Alm.append(Alm_grid_l3['Alm'])# Note: ftype = 'gauss' or 'gate' depends on the Alm grid file. No need to specify it	
	if lmax >=1:
		return Alm_grid_l1['theta'], Alm_grid_l1['delta'], Alm
	else:
		print('Warning: lmax<1. No grid can be retrieved')
		return [], [],[]

def numax_fct(Mass=1, Radius=1, Teff=5777.):
	numax_sun=3100.
	Teff_sun=5777.
	numax=numax_sun * Mass/(Radius**2 *np.sqrt(Teff/Teff_sun))
	return numax

def update_min(aj, aj_ref):
	if aj_ref > aj:
		return aj
	else:
		return aj_ref

def update_max(aj, aj_ref):
	if aj_ref < aj:
		return aj
	else:
		return aj_ref

#Assumptions: nu_max is derived from the requested Dnu_star parameter using the relation from Stello+2009. 
#	Dnu ~ 0.263*numax^0.77 (no uncertainty implemented here)
def numax_from_stello2009(Dnu_star):
	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	nu_max=np.power(10, np.log10(Dnu_star/beta0)/beta1)
	return nu_max

def Dnu_from_stello2009(numax):	
	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	Dnu=beta0*np.power(numax, beta1)
	return Dnu

def eval_aj_mean_range(numax=None, Dnu=None, Mass=None, Radius=None, Teff=None, epsilon_nl=5*1e-4, a1=1000, dir_core='/Users/obenomar/tmp/test_a2AR/', verbose=True):
	'''
	Evaluate the full range of possible of aj E [2,4] coefficients for a given star of Mass=Mass, Radius=Radius and Teff=Teff 
	and given an activity level epsilon_nl and a rotation a1
	If one wants a purely observational approach, it can provides directly numax (instead of M, R, Teff) and/or Dnu
	If only numax is provided, Dnu is determined using scaling relation with numax from Stello+2009. Same if only Dnu is provided
	This is usefull to set the range of simulations for aj and for setting priors on aj odds coefficient while fitting
	Note: This function **only** considers the mean value of aj over l and at numax
	'''
	R_sun=6.96342e5 #in km
	M_sun=1.98855e30 #in kg
	if numax != None and Dnu == None:
		if verbose == True:
			print('Calculating Dnu(numax)...')
		Dnu=Dnu_from_stello2009(numax)
		if verbose == True:
			print('Dnu = ', numax)
		eta0=eta0_fct(Dnu=Dnu)
	#
	if numax == None and Dnu != None:
		if verbose == True:
			print('Calculating numax(Dnu)...')
		numax=numax_from_stello2009(Dnu)
		eta0=eta0_fct(Dnu=Dnu)
	#
	if numax != None and Dnu != None:
		if verbose == True:
			print('Using numax and Dnu...')
		eta0=eta0_fct(Dnu=Dnu)
	#	
	if numax == None and Dnu == None and Mass != None and Radius != None and Teff !=None:
		if verbose == True:
			print('Using M,R,Teff...')
		numax=numax_fct(Mass=Mass, Radius=Radius, Teff=Teff)
		Volume=np.pi*(Radius*R_sun)**3 * 4./3
		rho=Mass*M_sun/Volume * 1e-12 # Conversion into g/cm-3
		#print('rho =', rho)
		eta0=eta0_fct(rho=rho)
	if numax == None and Dnu == None and (Mass == None or Radius == None or Teff == None):
		print('Error in eval_aj_mean_range: You need to provide (M, R, Teff) or numax and/or Dnu')
		print('                             Please read the description of the function before using it')
		exit()
	
	print('numax =', numax)
	#print('eta0  =', eta0)
	lmax=2 # We do not consider a6... we evaluate the range for a2 and a4 only
	dir_grids=dir_core+"/grids/gate/0.25deg_resol/" # High precision grid
	gridfiles=['grid_Alm_1.npz', 'grid_Alm_2.npz', 'grid_Alm_3.npz'] # Fine grids
	#
	theta, delta, Alm= load_Alm_grids(dir_grids, gridfiles, lmax=lmax)
	#
	a2_min=99999
	a2_max=-99999
	a4_min=99999
	a4_max=-99999
	a6_min=99999
	a6_max=-999999
	for j in range(len(theta)):
		for k in range(len(delta)):
			aj=np.zeros(6)
			for l in range(1,lmax+1):
				nu_nlm=[]
				for m in range(-l, l+1):
					# perturvation from AR		
					dnu_AR=numax*epsilon_nl*Alm[l-1][m+l][j,k] # All Alm(theta, delta) are in [j,k]
					# perturbation from CF
					dnu_CF=eta0*numax * (a1*1e-9)**2 * Qlm(l,m)
					nu_nlm.append(numax + dnu_AR + dnu_CF)
					#print('(j,k) = (', j, ',', k,')  theta =', theta[j]*180./np.pi, "  delta =", delta[k]*180./np.pi, " dnu_AR =", dnu_AR*1000)
				a=eval_acoefs(l, nu_nlm)
				# Averaging over l
				for o in range(len(aj)):
					if o < 2:
						aj[o]=aj[o] + a[o]/lmax
					if o >=2 and o < 4:
						if l>=2:
							aj[o]=aj[o] + a[o]/(lmax-1)
					if o >=3 and o <6:
						if l>=3:
							aj[o]=aj[o] + a[o]/(lmax-2)
#					if j==10 and k==10: # For debug only
#						print('a[{}] = {}'.format(o, a[o]))
#			if j==10 and k==10 and l==lmax: # For debug only
#				print('aj = ', aj)
#				exit()
			a2_min=update_min(aj[1], a2_min)
			a2_max=update_max(aj[1], a2_max)
			a4_min=update_min(aj[3], a4_min)
			a4_max=update_max(aj[3], a4_max)
			a6_min=update_min(aj[5], a6_min)
			a6_max=update_max(aj[5], a6_max)
	return [a2_min*1e3, a2_max*1e3], [a4_min*1e3, a4_max*1e3], [a6_min*1e3, a6_max*1e3] # aj converted into nHz


def aj_mean_range_v1():
	'''
	This function is to explore specific range for some stars. It can for example be used to evaluate 
	The prior that has to be used before fitting a star for which M, R, Teff is known.
	'''
	R_sun=6.96342e5 #in km
	M_sun=1.98855e30 #in kg
	V=np.pi*R_sun**3 * 4./3
	rho_sun=M_sun/V * 1e-12 # Conversion into g/cm-3
	# Table of M, R derived from Bellinger+2019 (https://www.aanda.org/articles/aa/pdf/2019/02/aa34461-18.pdf)
	# Only Sun, Min and Max of his range are kept
	# Teff are comming from Pinsonault+2012 or Davies+2016 (or Campante+2015)

	'''
	# USE THIS PART TO STUDY THE MAXIMUM POTENTIAL RANGE OF a2 and a4 FOR KEPLER OBSERVATIONS
	KIC=['6278762', 'Sun'    , '6679371' , '7199397' ]
	M=[0.753      ,   1      ,  1.570    , 1.308]
	R=[0.7588     ,   1      ,  2.222    , 2.526]
	rho=[1.741    , rho_sun  , 0.1432    ,  0.0821]
	T=[5046       , 5777.    , 6551      ,  5821.]
	epsilon_nl=5e-4
	#a1=[500,  1000, 1500, 2000, 2500]
	a1=[420]
	'''
	# USE THIS PART FOR GETTING SPECIFIC EXPECTED MAXIMUM RANGE FOR THE REFERENCE STAR OF SIMULATION: 16 CYG A
	# Data from Bellinger+2019 for M, R, rho. Teff from Ramirez+2009 (and Metcalfe+2012)
	# a1 is set as per the simulations and epsilon_nl is the solar value
	
	KIC=['12069424']
	M=[1.056]
	R=[1.213]
	T=[5825]
	rho=[1.531]
	a1=[450]
	epsilon_nl=5e-4 # Solar Value
	print('epsilon_nl :', epsilon_nl)
	
	''''
	KIC=['12069449']
	M=[1.07]
	R=[1.127]
	T=[5750]
	#rho=[1.531]
	a1=[500]
	epsilon_nl=5e-4 # Solar Value
	'''
	#for j in range(len(M)):
	for j in range(0,1):
		print('j=', j)
		print(colored('KIC =', 'red'), KIC[j])
		print('Mass: ', M[j],  '   Radius: ', R[j],   '    Teff: ', T[j])
		for dnu in a1:
			print(colored('a1         :'+str(dnu), 'red'))
			a2_range, a4_range, a6_range=eval_aj_mean_range(Mass=M[j], Radius=R[j], Teff=T[j], epsilon_nl=epsilon_nl, a1=dnu)
			print(colored('a2_range:', 'red'), a2_range, ' (nHz) | ', colored(100*np.abs(a2_range[0])/dnu, 'blue'),  ' , ', colored(np.abs(100*a2_range[1])/dnu, 'blue'), '  (%)')
			print(colored('a4_range:', 'red'), a4_range, ' (nHz) | ', colored(100*np.abs(a4_range[0])/dnu, 'blue'),  ' , ', colored(np.abs(100*a4_range[1])/dnu, 'blue'), '  (%)')
			print(colored('a6_range:', 'red'), a6_range, ' (nHz) | ', colored(100*np.abs(a6_range[0])/dnu, 'blue'),  ' , ', colored(np.abs(100*a6_range[1])/dnu, 'blue'), '  (%)')
			print(colored('----', 'red'))


def aj_mean_range_v2(epsilon_nl=5e-4, IDs_file='../../data/id_list.txt', 
					 Bellinger_tabA1='../../External_data/Legacy_params_Bellinger2019.txt',
					 Pinsonnault_tab7='../../External_data/Pinsonnault2012_table7.dat',
					 Benomar_tab='../../External_data/Level2-Outputs/a1a3_summary.txt',
					 output_file='../../data/aj_range_RESULTS.txt',
					 dir_core =None):
	'''
	This function is to explore specific range for some stars. It can for example be used to evaluate 
	The prior that has to be used before fitting a star for which M, R, Teff is known.
	It reads a file that contain ID, M, R, T, rho and a1 for each star and then returns in the screen a full list of the ranges for aj
	'''
	#dir_core='/Users/obenomar/Work/tmp/paper2_asphericity/'
	if dir_core == None:
		dir_core=os.cwd()

	R_sun=6.96342e5 #in km
	M_sun=1.98855e30 #in kg
	V=np.pi*R_sun**3 * 4./3
	rho_sun=M_sun/V * 1e-12 # Conversion into g/cm-3
	# Reading the IDs within the provided file
	print('Reading the IDs within the provied file: ', IDs_file)
	myID_list=read_id_list(IDs_file)
	# Table of M, R derived from Bellinger+2019 (https://www.aanda.org/articles/aa/pdf/2019/02/aa34461-18.pdf)
	print('Reading Bellinger et al. 2019 table A1...')
	data_bellinger, header_bellinger=read_matrix_tab(Bellinger_tabA1)
	i_M=2
	i_R=3
	# Teff are comming from Pinsonault+2012
	print('Reading Pinsonnault et al. 2012 table 7...')
	data_pinsonnault, header_pinsonnault=read_matrix_tab(Pinsonnault_tab7)
	i_Teff=1
	# Benomar et al 2018 for a1 table
	print('Reading Benomar et al. 2018 summary table of a1 and a3...')
	data_benomar, header_benomar=read_matrix_tab(Benomar_tab)
	i_a1=1
	# Cross referencing
	print('Attempting to match ID of the published tables with the one provided by the user...')
	print('   -  Verifying that all ID from the user exist...')
	for myID in myID_list:
		id_pos, label, err=match_id_within_lists(myID, data_bellinger, data_pinsonnault, data_benomar, verbose=True)
	if err == True:
		raise ValueError(" Error while trying to match star catalogs")
	print('    - Calculating the range of a2 and a4...')
	txtforfile='#Ranges defined using acoefs2activity/eval_aj_range.py\n'
	txtforfile=txtforfile + '#Ranges are provided in microHz\n'
	txtforfile=txtforfile + '#epsilon_nl = {}\n'.format(epsilon_nl)
	txtforfile=txtforfile + '#ID             a2_min           a2_max          a4_min         a4_max\n'
	for j in range(len(myID_list)):
		print(colored('[{}/{}] ID ='.format(j+1, len(myID_list)), 'red'), colored(myID_list[j], 'red'))
		id_pos, label=match_id_within_lists(myID_list[j], data_bellinger, data_pinsonnault, data_benomar, verbose=False)
		M=data_bellinger[id_pos[0], i_M]
		R=data_bellinger[id_pos[0],i_R]
		Teff=data_pinsonnault[id_pos[1], i_Teff]
		a1=data_benomar[id_pos[2], i_a1]
		print('    Mass: ', M,  '   Radius: ', R,   '    Teff: ', Teff, '     a1: ', a1)
		a2_range, a4_range, a6_range=eval_aj_mean_range(dir_core=dir_core, Mass=M[0], Radius=R[0], Teff=Teff[0], epsilon_nl=epsilon_nl, a1=a1[0], verbose=False)
		print(colored('a2_range:', 'red'), a2_range, ' (nHz) | ', colored(100*np.abs(a2_range[0])/a1[0], 'blue'),  ' , ', colored(np.abs(100*a2_range[1])/a1[0], 'blue'), '  (%)')
		print(colored('a4_range:', 'red'), a4_range, ' (nHz) | ', colored(100*np.abs(a4_range[0])/a1[0], 'blue'),  ' , ', colored(np.abs(100*a4_range[1])/a1[0], 'blue'), '  (%)')
		#print(colored('a6_range:', 'red'), a6_range, ' (nHz) | ', colored(100*np.abs(a6_range[0])/a1[0], 'blue'),  ' , ', colored(np.abs(100*a6_range[1])/a1[0], 'blue'), '  (%)')
		print(colored('----', 'red'))
		#txtforfile=txtforfile + '{0}  {1:14.5f}  {2:14.5f}  {3:14.5f}  {4:14.5f}  {5:14.5f}  {6:14.5f}\n'.format(myID_list[j], a2_range[0]/1e3, a2_range[1]/1e3, a4_range[0]/1e3, a4_range[1]/1e3, a6_range[0]/1e3, a6_range[1]/1e3)
		txtforfile=txtforfile + '{0}  {1:14.5f}  {2:14.5f}  {3:14.5f}  {4:14.5f}\n'.format(myID_list[j], a2_range[0]/1e3, a2_range[1]/1e3, a4_range[0]/1e3, a4_range[1]/1e3)
	f=open(output_file, 'w')
	f.write(txtforfile)
	f.close()

def aj_mean_range_v3(epsilon_nl=5e-4, IDs_file='../../data/id_list_Kamiaka2018_NoLegacy.txt', 
					 dir_tamcmc_data="../../data/REFERENCE_RUN_1001/Kamiaka_stars/14Aug2023/outputs/",
					 output_file='../../data/aj_range_Kamiaka2018_NoLegacy.txt',
					 dir_grids="../../", cpp_path="../../../../TAMCMC_bin_1.85.1_AMD/"):
	'''
	This function is to explore specific range for some stars. It can for example be used to evaluate 
	The prior that has to be used before fitting a star for which M, R, Teff is known.
	It returns in the screen a full list of the ranges for aj
	Contrary to the v2, a1 is read from a previous analysis of a star with (preferably only) a1 as a aj parameter.
	It also use that previous analysis to measure Dnu and estimate numax.
	This is usefull if existing table of M,R and Teff do not necessarily exist. But may be less precise while more data driven
	'''

	txtforfile='#Ranges defined using acoefs2activity/eval_aj_range.py with aj_mean_range_v3\n'
	txtforfile=txtforfile + '#Ranges are provided in microHz\n'
	txtforfile=txtforfile + '#epsilon_nl = {}\n'.format(epsilon_nl)
	txtforfile=txtforfile + '#ID             a2_min           a2_max          a4_min         a4_max\n'
	
	R_sun=6.96342e5 #in km
	M_sun=1.98855e30 #in kg
	V=np.pi*R_sun**3 * 4./3
	rho_sun=M_sun/V * 1e-12 # Conversion into g/cm-3
	# Reading the IDs within the provided file
	print('Reading the IDs within the provied file: ', IDs_file)
	myID_list=read_id_list(IDs_file)
	splitting_labels=["a1", "a1_0", "splitting_a1"] # accepted keywords for the splitting
	ID_err=[]
	j=0
	splitting_labels=["a1", "a1_0", "splitting_a1", "sqrt(splitting_a1).cosi", "sqrt(splitting_a1).sini"]
	for ID in myID_list:
		print(colored('[{}/{}] ID ='.format(j+1, len(myID_list)), 'red'), colored(myID_list[j], 'red'))
		process_name=get_dirs_list(dir_tamcmc_data, ID, subdirectory_only=True) # Look for directory which has the name of the ID within it
		if process_name == []:
			print("Could not find the star within the directory ", dir_tamcmc_data)
			print("   Star ID: ", ID)
			raise ValueError("Error: Cannot find the data associated to the provided ID")
		elif len(process_name) != 1:
			print("Found multiple stars with the same name")
			raise ValueError("Error: We cannot handle multiple similar names")
		else:
			process_name=process_name[0]
		print('Reading samples previous MCMC analysis')
		samples, labels, isfixed, plength=bin2txt(dir_tamcmc_data, process_name, phase='A', chain=0, 
			first_index=0, last_index=-1, period=20, single_param_index=-1,
			erase_tmp=True, cpp_path=cpp_path, cpp_version="1.85.0", outdir='../../tmp/', get_plength=True)
		inds_splitting=[]
		for seek_label in splitting_labels:
			for k in range(len(labels)):
				if labels[k] == seek_label:
					inds_splitting.append(k)
		print('Determining a1...')
		if inds_splitting == []:
			print(" Error: inds_splitting is empty !")
			raise ValueError("Cannot proceed without a1")
		if len(inds_splitting) > 2:
			print(" Error : inds_splitting contains multiple possible indexes")
			print("       labels[inds_splitting] : ")
			for ind in inds_splitting:
				print("       ", labels[ind])
			print("Cannot proceed")
			raise ValueError("Cannot proceed due to multiple splitting values")
		if len(inds_splitting) == 2:
			if ("sqrt(splitting_a1).cosi" in labels[inds_splitting[0]] and "sqrt(splitting_a1).sini" in labels[inds_splitting[1]]) or \
				("sqrt(splitting_a1).cosi" in labels[inds_splitting[1]] and "sqrt(splitting_a1).sini" in labels[inds_splitting[0]]):
				print(" Conversion from 'sqrt(splitting_a1).cosi', 'sqrt(splitting_a1).sini' required...")
				a1_samples=samples[:,inds_splitting[0]]**2 + samples[:,inds_splitting[1]]**2
				a1, errors_a1, stats_a1=uncertainties(a1_samples)
			else:
				raise ValueError("Error: Conversion of a1 impossible for an unknown reason. Debug required")
		if len(inds_splitting) == 1:
			a1, errors_a1, stats_a1=uncertainties(samples[:,inds_splitting[0]])
		a1=a1*1000 # convert to nHz
		errors_a1=errors_a1*1000
		stats_a1=stats_a1*1000
		a3_range=[-0.2*a1, 0.2*a1] # a3 is max 20% of a1
		print(" a1 = {}    err+ = {}     err- = {}".format(a1, errors_a1[1],errors_a1[0]))
		print('Determining Dnu...')
		nu_nl=get_nu_samples(samples, plength, verbose=False, step_factor=1, first_sample=0)
		Dnu_samples=get_Dnu_samples(nu_nl[0,:,:])
		Dnu, errors_Dnu, stats_Dnu=uncertainties(Dnu_samples)
		print(" Dnu = {}    err+ = {}     err- = {}".format(Dnu, errors_Dnu[1],errors_Dnu[0]))
		print('    - Calculating the range of a2 and a4...')
		print('    Dnu: ', Dnu  ,'     a1: ', a1)
		a2_range, a4_range, a6_range=eval_aj_mean_range(dir_core=dir_grids, Dnu=Dnu, epsilon_nl=epsilon_nl, a1=a1, verbose=False)
		print(colored('a2_range:', 'red'), a2_range, ' (nHz) | ', colored(100*np.abs(a2_range[0])/a1, 'blue'),  ' , ', colored(np.abs(100*a2_range[1])/a1, 'blue'), '  (%)')
		print(colored('a3_range:', 'red'), a3_range, ' (nHz) | ', colored(100*np.abs(a3_range[0])/a1, 'blue'),  ' , ', colored(np.abs(100*a3_range[1])/a1, 'blue'), '  (%)')
		print(colored('a4_range:', 'red'), a4_range, ' (nHz) | ', colored(100*np.abs(a4_range[0])/a1, 'blue'),  ' , ', colored(np.abs(100*a4_range[1])/a1, 'blue'), '  (%)')
		#print(colored('a6_range:', 'red'), a6_range, ' (nHz) | ', colored(100*np.abs(a6_range[0])/a1, 'blue'),  ' , ', colored(np.abs(100*a6_range[1])/a1, 'blue'), '  (%)')
		print(colored('----', 'red'))
		txtforfile=txtforfile + '{0}  {1:14.5f}  {2:14.5f}  {3:14.5f}  {4:14.5f}  {5:14.5f}  {6:14.5f}\n'.format(myID_list[j], a2_range[0]/1e3, a2_range[1]/1e3,  a3_range[0]/1e3, a3_range[1]/1e3, a4_range[0]/1e3, a4_range[1]/1e3)	
		j=j+1	
	f=open(output_file, 'w')
	f.write(txtforfile)
	f.close()
	if ID_err !=[]:
		print("Failed stars : ")
		for ID in ID_err:
			print("     ", ID)

def aj_range4modelfiles(range_file='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/aj_range_RESULTS.txt'):
	'''
		A small routine that read a file created by aj_mean_range_v2() and 
		show it is in a handy way for copy/past inside a .model file
	'''
	data,header=read_matrix_tab(range_file)
	for i in range(len(data[:,0])):
		print(colored('KIC ' + str(int(data[i,0])),'red'))
		print(colored('------------', 'red')) 
		a2_line_0=' a2_0                     Uniform          {0:6.4f}       {1:6.4f}       {2:6.4f}'.format(0.0000, data[i,1], data[i,2])
		a2_line_1=' a2_1                     Fix              0.000 '
		a4_line_0=' a4_0                     Uniform          {0:6.4f}       {1:6.4f}       {2:6.4f}'.format(0.0000, data[i,3], data[i,4])
		a4_line_1=' a4_1                     Fix              0.000 '
		print(a2_line_0)
		print(a2_line_1)
		print(' ')
		print(a4_line_0)
		print(a4_line_1)
		print(colored('-------- ', 'red')) 
		key=input('Press any key to proceed to the next star...')
 

##aj_mean_range_v1()
##aj_mean_range_v2(epsilon_nl=1e-3)
#aj_range4modelfiles()