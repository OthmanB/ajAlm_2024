import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
from  scipy.io import readsav
#from plotly.subplots import make_subplots
from matplotlib import gridspec
import matplotlib.patches as mpatches
from scipy import interpolate
import matplotlib.cm as cm
import seaborn as sns # Very good for statistical plotting
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from read_aj_mcmcdata import *

def get_a1ovGamma_atnumax(loga1_ov_gamma_file, nu_l0, numax, log_results=False, fast=False):
	'''
		This program gets as input 
			- A sav file that contains the samples LOGA1_OV_WIDTH (generated with my IDL postMCMC program)
		for each of modes of degree l and for all n
			- The list of frequencies nu_l0=nu(n,l=0, m=0) that may be extracted from the IDL postMCMC files
		or from the model file 
			- A value of numax
		From this, and using l=0 values, it determines the samples of a1_ov_width_at_numax by function interpolation 
		log_results: If True, return the log of the ratio a1ovGamma instead of the ratio itself
		fast : reduce by 'ratio' the number of samples (the first 10%) to reduce the computation time
	'''
	d=readsav(loga1_ov_gamma_file)
	Nd=len(d['loga1_ov_width'][0,:,0]) # (samples, n, l)
	if fast == False:
		Nsamples=len(d['loga1_ov_width'][:,0,0])
	else:
		ratio=20
		Nsamples=int(len(d['loga1_ov_width'][:,0,0])/ratio)
		cpt=0
	#xfit=np.linspace(0, Nd-1, Nd)
	xfit=nu_l0
	loga1_ov_gamma_at_numax=np.zeros(Nsamples)
	for i in range(Nsamples):
		if fast == False:
			interpol=interpolate.interp1d(xfit, d['loga1_ov_width'][i,:,0], kind='cubic')
		else: # In fast mode, we take one sample every 10
			interpol=interpolate.interp1d(xfit, d['loga1_ov_width'][cpt,:,0], kind='cubic')
			cpt=cpt+ratio
		loga1_ov_gamma_at_numax[i]=interpol(numax)
	if log_results == False:
		return np.exp(loga1_ov_gamma_at_numax)
	else:
		return loga1_ov_gamma_at_numax


def bias_pdf(samples, xlabels, text, xtrue=None, Ncols=3, Nrows=3, file_out='bias'):
	'''
		Function that receives a 2D array containing samples ('samples') for each parameters
		For each parameter, a pdf is shown in a single-page that contains (Ncols,Nrows) fields.
		xlabel: Defines the name of the x-axis (name of the parameter)
		ylabel: Defines the name of the y-axis (PDF)
		file_out: root name/directory for the image
		text: Any kind of text that we would like to be embeded inside the image of each field
		xtrue: If set, provides the value of the true input. It must be of same size as the number of parameters
	'''
	bins=30#75
	Nparams=len(samples[:,0]) # Ninc
	#	
	# Evaluate uncertainties using the samples
	errors=np.zeros((2,Nparams))
	med=np.zeros(Nparams)
	for i in range(Nparams):
		stats = np.percentile(samples[i,:], [16, 50, 84])
		#print(' stats[',i,'] = ', stats)
		errors[0,i]=stats[1] - stats[0]
		errors[1,i]=stats[2] - stats[1]
		med[i]=stats[1]
	#
	cpt=0
	fig_1d, ax = plt.subplots(Ncols,Nrows, figsize=(12, 6))
	for i in range(Ncols):
		for j in range(Nrows):
			h=ax[i,j].hist(samples[cpt,:],bins=bins)
			#pos=np.where(np.bitwise_and(h[1] >= med[cpt] - errors[0,cpt], h[1] <= med[cpt] + errors[1,cpt]))
			#ax[i,j].fill(h[1][pos], h[0][pos], color='gray')
			if j == Nrows: # Put the label only once, at the bottom of the plots
				ax[i,j].set_xlabel(xlabels)
			if i == Ncols: # Put the label only on the left side of the plots
				ax[i,j].set_ylabel("PDF")
			#print("[{}]   med[{}] = {}".format(cpt, cpt, med[cpt]))
			#print("     xtrue[{}] = {}".format(cpt, xtrue[cpt]))		
			ax[i,j].axvline(x=med[cpt], color='blue', linestyle='-')
			ax[i,j].axvline(x=med[cpt]-errors[0,cpt], color='gray', linestyle='--')
			ax[i,j].axvline(x=med[cpt]+errors[1,cpt], color='gray', linestyle='--')
			ax[i,j].axvline(x=xtrue[cpt], color='red')
			ax[i,j].text(0.8,0.8, text[cpt], fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax[i,j].transAxes)
			cpt=cpt+1
	fig_1d.savefig(file_out+'_pdf.jpg')
	plt.close('all')

def read_combi(file, filter_HNR=None):
	''' 
		Program dedicated to the reading of the Combination files that are generated by the Spectra-Simulator program
	'''
	pnames_done=False
	f=open(file, 'r')
	txt=f.read()
	f.close()
	txt=txt.split('\n')
	data=[]
	for t in txt:
		s=t.split()
		if s != '' and s !=[]:
			if s[0] != '':
				if s[0] == 'model_name=':
					model_name=s[1]
				if s[0] == '#': # This is the labels
					if pnames_done == False:
						params_names=s[1:]#.split()
						pnames_done=True
					else:
						print('Error: For some reason, there was more than one indicator of label. The label indicator is a #')
						print('       Be sure not to edit the Combinations.txt generated by the Spectra-Simulator program')
						exit()
				else: # This must be the data
					if pnames_done == True:
						data.append(s)
	if filter_HNR != None:
		HNRs=lookup_combi('HNR', data, params_names, return_index=False)
		posOK=np.where(np.asarray(HNRs, dtype=float) == filter_HNR)
		if len(posOK[0]) == 0:
			print(' Error: filter_HNR = {} but HNRs do not contain such a value'.format(filter_HNR))
			print('        Debug required: ')
			print('        HNRs = {}'.format(HNRs))
			exit()
		datanew=[]
		for ind in range(len(posOK[0])):
			datanew.append(data[posOK[0][ind]][:])
		return datanew, params_names, model_name			
	else:
		return data, params_names, model_name

def lookup_combi(key, combi_data, param_names, return_index=False):
	'''
		Look for a keyword 'key' inside param_names and returns the column of data associated to that key
		If the key does not exist, an error is returned.
		return_index: If True, returns only the index at which the key exists instead of the column
	'''
	index=param_names.index(key)
	tab=[]
	if return_index == False:
		for c in combi_data:
			tab.append(c[index])
		return tab
	else:
		return index


def bias_analysis_v1(MCMCdir, combi_file, fileout='plot'):
	'''
		Note on v1 version: This version shows a single Scenario (eg Polar case OR Equatorial case). In that sense, it 
							can only accept one MCMCdir and one combi_file
							To Show multiple scenarii in a single plot, use the bias_analysis_v2()
		This program gather informations on aj from a rootdirectory that contains an ensemble of simulations and 
		compare them with true inputs provided by a 'Combinations.txt' file. 
		The Combination file is generated by the Spectrum-Simulator as per specified by the configuration of the simulator:
		(see https://github.com/OthmanB/Spectra-Simulator-C)
		The outputs are:
			- A table of all the aj compared to the true inputs [NOT IMPLEMENTED]
			- Plots showing aj in function of inclination
	'''
	confidence=[2.25,16,50,84,97.75]
	jmax=6
	unit_nHz=1000 # Multiplicator to get units in nHz
	#
	print('A. Reading the Combinations file...')
	combi_data, param_names, model_name=read_combi(combi_file)
	Ncombi=len(combi_data)
	#print('B. Identify the directories in the rootdir...')
	#list_dirs=get_dirs_list(MCMCdir) 
	#print(list_dirs)
	print('B. Extracting aj and inclinations for each directory of that has a name that match the ID name in the combi_file... ')
	
	aj_data=np.zeros((Ncombi,jmax, len(confidence)))
	inc_data=np.zeros((Ncombi, len(confidence)))
	for i in range(Ncombi):
		stardir=MCMCdir + '/' + combi_data[i][0] + '.0/' # That should be a sequence of numbers that identify the simulation
		print('  --------')
		print('  star:', combi_data[i][0])
		print('  --------')
		aj_stats, inc_stats, aj_spl, inc_spl, Nsamples=get_aj_inc_star(stardir, confidence=confidence)
		aj_data[i,:,:]=aj_stats*unit_nHz
		inc_data[i,:]=inc_stats
		if i == 0:
			aj_samples=np.zeros((Ncombi, jmax, Nsamples))
			inc_samples=np.zeros((Ncombi, Nsamples))
		aj_samples[i,:,:]=aj_spl*unit_nHz
		inc_samples[i,:]=inc_spl
	print('C. Plots of bias...')
	inc_true=np.asarray(lookup_combi('i', combi_data, param_names), dtype=float)
	inc0=30 # Reject any plot of bias below inc0
	posOK=np.where(inc_true >= inc0)
	Ninc=len(inc_true)
	for j in range(0,jmax):
		print(" j=", j+1, ":")
		print("    1. Bias as a function of inclination")
		yr=[-1000, 1000] # Default yrange
		if j+1 != 1:
			aj_true=np.asarray(lookup_combi('a'+str(j+1), combi_data, param_names), dtype=float)*unit_nHz
			if j+1 == 2:
				yr=[-250, 250]
			if j+1 == 4:
				yr=[-20, 20]				
		else: # a1_true is not an explicit parameter. There is two choice: Use a1ovGamma*Gamma_at_numax or set it manually.
			a1ovGamma=np.asarray(lookup_combi('a1ovGamma', combi_data, param_names), dtype=float)
			Gamma_at_numax=np.asarray(lookup_combi('Gamma_at_numax', combi_data, param_names), dtype=float)
			a1_true=a1ovGamma * Gamma_at_numax*unit_nHz
			aj_true=np.zeros(Ninc) + a1_true
			yr=[np.min(aj_true[posOK])*0.1, np.max(aj_true[posOK])*1.1]
		minfill=min(yr)
		fig ,ax= plt.subplots(tight_layout=True)
		ax.set_xlabel('Inclination (deg)')
		#ax.set_ylabel('$a_'+str(j+1)+' - a_'+str(j+1) + '^{(true)}$ (nHz)')
		ax.set_ylabel('$a_'+str(j+1) + '$ (nHz)')
		ax.plot(inc_true, aj_true, linestyle='dashed',dashes=(5, 10)) # Horizontal line
		#yr=[min([np.min(aj_true[posOK]), np.min(aj_data[posOK,j,2])]), max([np.max(aj_true[posOK]), np.max(aj_data[posOK,j,2])])]
		ax.set_ylim(yr[0], yr[1])
		ax.set_xlim(-5, 95)
		yerr=make_error_from_stats(aj_data[:,j,:])
		xerr=make_error_from_stats(inc_data[:,:])
		#ax.errorbar(inc_data[posOK,2], aj_data[posOK,j, 2], xerr=xerr[posOK], yerr=yerr[posOK])
		for i in range(Ninc):
			if inc_true[i] >= inc0:
				ax.annotate("", xy=(inc_data[i,2], aj_data[i,j,2]), xytext=(inc_true[i], aj_true[i]), va="center", ha="center", arrowprops={"arrowstyle": "-|>"})
		plt.fill_between([0, inc0], [yr[1],yr[1]], minfill, color='gray')
		plt.savefig(fileout+'_a'+str(j+1)+'.png', dpi=300)
		plt.close('all')
		#
		print("    2. Bias for each aj (PDF)")
		xlabel='$a_'+str(j+1) + '$ (nHz)'
		i0=1 # Ignore the inclination=0 to have a (3,3) = 9 Plot
		text=[]
		for i in range(i0,len(inc_true)):
			text.append('inc='+str(inc_true[i]))
		bias_pdf(aj_samples[i0:,j,:], xlabel, text, xtrue=aj_true[i0:], Ncols=3, Nrows=3, file_out=fileout + '_bias_a'+str(j+1))

def bias_analysis_v2(MCMCdir, combi_files, labels=None, fileout='plot', col_arrows=None):
	'''
		Note on v2 version: This version shows multiple scenarii in a single plot, use the bias_analysis_v2()
							To Show a single Scenario (eg Polar case OR Equatorial case) use bias_analysis_v1()
		This program gather informations on aj from a rootdirectory that contains an ensemble of simulations and 
		compare them with true inputs provided by a 'Combinations.txt' file. 
		The Combination file is generated by the Spectrum-Simulator as per specified by the configuration of the simulator:
		(see https://github.com/OthmanB/Spectra-Simulator-C)
		The outputs are:
			- A table of all the aj compared to the true inputs [NOT IMPLEMENTED]
			- Plots showing aj in function of inclination
			- Plots showing aj PDFS with an indicator of the true input
	'''
	confidence=[2.25,16,50,84,97.75]
	jmax=6
	unit_nHz=1000 # Multiplicator to get units in nHz
	inc0=30 # Reject any plot of bias below inc0
	#
	print('A. Reading the Combinations files...')
	combi_data=[]
	param_names=[]
	model_name=[]
	Ncombi=[]
	for c in combi_files:
		c_data, p_names, m_name=read_combi(c)
		combi_data.append(c_data)
		param_names.append(p_names)
		model_name.append(m_name)
		Ncombi.append(len(c_data))
	Nscenario=len(Ncombi)
	if col_arrows == None:
		col_arrows=np.repeat('Black', Nscenario)	
	if labels != None:
		patches=[]
		for k in range(0, Nscenario):
			patches.append(mpatches.Patch(color=col_arrows[k], label=labels[k]))

	col_l=['black', 'dimgray']
	col_lines=[['black', 'black']]
	for j in range(1, jmax):
		col_lines.append(col_l)
	print('B. Extracting aj and inclinations for each directory of that has a name that match the ID name in the combi_files... ')
	aj_data=[]
	inc_data=[]
	aj_samples=[]
	inc_samples=[]
	for k in range(Nscenario):
		print('  ********** SCENARIO {} **********  '.format(k))
		print('  MCMCdir =', MCMCdir[k])
		print('  *********************************  ')
		aj_data.append(np.zeros((Ncombi[k],jmax, len(confidence))))
		inc_data.append(np.zeros((Ncombi[k], len(confidence))))
		for i in range(Ncombi[k]):
			stardir=MCMCdir[k] + '/' + combi_data[k][i][0] + '.0/' # That should be a sequence of numbers that identify the simulation
			print('  --------')
			print('  star:', combi_data[k][i][0])
			print('  --------')
			aj_stats, inc_stats, aj_spl, inc_spl, Nsamples =get_aj_inc_star(stardir, confidence=confidence)
			aj_data[k][i,:,:]=aj_stats*unit_nHz
			inc_data[k][i,:]=inc_stats
			if i == 0:
				aj_samples.append(np.zeros((Ncombi[k], jmax, Nsamples)))
				inc_samples.append(np.zeros((Ncombi[k], Nsamples)))
			aj_samples[k][i,:,:]=aj_spl*unit_nHz
			inc_samples[k][i,:]=inc_spl

	print('C. Plots of bias...')
	inc_true=np.asarray(lookup_combi('i', combi_data[k], param_names[k]), dtype=float)
	posOK=np.where(inc_true >= inc0)
	Ninc=len(inc_true)
	for j in range(0,jmax):
		print(" j=", j+1, ":")
		yr=[-1000, 1000] # Default yrange
		if j+1 == 1:
			yr=[900, 1100]
		if j+1 == 2:
			yr=[-250, 250]
		if j+1 == 4:
			yr=[-20, 20]		
		fig ,ax= plt.subplots(tight_layout=True)
		ax.set_xlabel('Inclination (deg)')
		ax.set_ylabel('$a_'+str(j+1) + '$ (nHz)')
		ax.set_ylim(yr[0], yr[1])
		ax.set_xlim(-1, 91)
		#fig.set_label(labels[k])
		#plt.legend(labels=labels, color=col_arrows, loc="upper left")
		for k in range(Nscenario):
			if j+1 != 1:
				aj_true=np.asarray(lookup_combi('a'+str(j+1), combi_data[k], param_names[k]), dtype=float)*unit_nHz		
			else: # a1_true is not an explicit parameter. There is two choice: Use a1ovGamma*Gamma_at_numax or set it manually.
				a1ovGamma=np.asarray(lookup_combi('a1ovGamma', combi_data[k], param_names[k]), dtype=float)
				Gamma_at_numax=np.asarray(lookup_combi('Gamma_at_numax', combi_data[k], param_names[k]), dtype=float)
				a1_true=a1ovGamma * Gamma_at_numax*unit_nHz
				aj_true=np.zeros(Ninc) + a1_true
			minfill=min(yr)
			yerr=make_error_from_stats(aj_data[k][:,j,:])
			xerr=make_error_from_stats(inc_data[k][:,:])
			#ax.errorbar(inc_data[posOK,2], aj_data[posOK,j, 2], xerr=xerr[posOK], yerr=yerr[posOK])
			ax.plot(inc_true, aj_true, linestyle='dashed',dashes=(5, 10), color=col_lines[j][k]) # Horizontal line
			for i in range(Ninc):
				if inc_true[i] >= inc0:
					ax.annotate("", xy=(inc_data[k][i,2], aj_data[k][i,j,2]), xytext=(inc_true[i], aj_true[i]), color=col_arrows[k], va="center", ha="center", arrowprops={"arrowstyle": "-|>", "color":str(col_arrows[k])})
			plt.fill_between([0, inc0], [yr[1],yr[1]], minfill, color='silver')
		if labels != None:
			ax.legend(handles=patches, loc="upper left")
		plt.savefig(fileout+'_a'+str(j+1)+'.png', dpi=300)
		plt.close('all')


def saturate_bias(bias, incs, inc0, extend_limit=0):
	'''
		Modify the bias vector in order to make it saturate to a 'limit' value specified
		by the maximum values of bias with the condition incs>inc0.
		By default, we saturate at the maximum value of the bias with incs>inc0 (extend_limit=0)
		But the user can request a higher value proportional to it, using extend_limit>0 (in fraction of limit).
	'''
	b=bias
	posOK = np.where(incs>=inc0)
	limit=np.abs(b[posOK]).max() # The limit is on the absolute value of the bias
	pos_p=np.where(b > limit) # Saturate on positive biases
	b[pos_p]=limit + limit*extend_limit
	pos_m=np.where(b < -limit) # Saturate on negative biases
	b[pos_m]=-limit - limit*extend_limit
	return b

def get_true_vals(j, combi_data, param_names, unit_nHz):
	'''
		Function that returns aj_true, Gamma_at_numax_true, a1ovGamma_true
		provided a table from combi_data
	'''
	a1ovGamma_true=np.asarray(lookup_combi('a1ovGamma', combi_data, param_names), dtype=float)
	Gamma_at_numax_true=np.asarray(lookup_combi('Gamma_at_numax', combi_data, param_names), dtype=float)
	Ntrue=len(a1ovGamma_true)
	if j != 1:
		aj_true=np.asarray(lookup_combi('a'+str(j), combi_data, param_names), dtype=float)*unit_nHz		
	else: # a1_true is not an explicit parameter. There is two choice: Use a1ovGamma*Gamma_at_numax or set it manually.
		a1_true=a1ovGamma_true * Gamma_at_numax_true*unit_nHz
		aj_true=np.zeros(Ntrue) + a1_true
	return aj_true, Gamma_at_numax_true, a1ovGamma_true

def bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout='bias', abs_err=False, saturate_colors=False, filter_HNR=None):
	'''
		Note on v3 version: This version creates a bias map using arrows in the (a1/Gamma, inc) space and 
							color/symbol code the information of the bias on aj (j=[1,6]) at point in the(a1/Gamma,inc) space
							It encodes a lot of information so it may be oversaturate
		This program gather informations on aj, a1, Gamma and inc from a rootdirectory that contains an ensemble of simulations and 
		compare them with true inputs provided by a 'Combinations.txt' file. 
		The Combination file is generated by the Spectrum-Simulator as per specified by the configuration of the simulator:
		(see https://github.com/OthmanB/Spectra-Simulator-C)
		Inputs:
			MCMCdir: A list of main directory containing all of the data for all the scenario considered
			combi_files: A list pointing towards all of the Combinations.txt files that define each considered scenario
			numax_star: numax of the star used for the simulation. Check the reference star model to know this
			fileout: rootnames of any file that is created on the process
			load_npz: If True, attempts to load a summary of all of the data, skipping the pre-processing phase. If False, it will create/overwrite the npz file
		The outputs are:
			- A table of all the aj compared to the true inputs [NOT IMPLEMENTED]
			- Plots showing aj in function of the ratio a1/Gamma and of inclination
	'''	
	confidence=[2.25,16,50,84,97.75]
	jmax=6
	unit_nHz=1000 # Multiplicator to get units in nHz
	inc0=35 # Specific treatment for known highly biased values when inc<inc0
	colorBAD='darkgray' # Color used for bias(inc<inc0) 
	#
	print('A. Reading the Combinations files...')
	combi_data=[]
	param_names=[]
	model_name=[]
	Ncombi=[]
	for c in combi_files:
		c_data, p_names, m_name=read_combi(c, filter_HNR=filter_HNR)
		combi_data.append(c_data)
		param_names.append(p_names)
		model_name.append(m_name)
		Ncombi.append(len(c_data))
	Nscenario=len(Ncombi)

	print('B. Extracting l=0 Frequencies, a1ovGamma, aj and inclinations for each directory of that has a name that match the ID name in the combi_files... ')
	aj_data=[]
	inc_data=[]
	a1ovGamma_data=[]
	aj_samples=[]
	inc_samples=[]
	a1ovGamma_samples=[]
	for k in range(Nscenario):
		print('  ********** SCENARIO {} **********  '.format(k))
		print('  MCMCdir =', MCMCdir[k])
		print('  *********************************  ')
		aj_data.append(np.zeros((Ncombi[k],jmax, len(confidence))))
		inc_data.append(np.zeros((Ncombi[k], len(confidence))))
		a1ovGamma_data.append(np.zeros(Ncombi[k]))
		for i in range(Ncombi[k]):
			stardir=MCMCdir[k] + '/' + combi_data[k][i][0] + '.0/' # That should be a sequence of numbers that identify the simulation
			print('  --------')
			print('  star:', combi_data[k][i][0])
			print('  --------')
			print('		(a) Extracting aj and inclinations...')
			aj_stats, inc_stats, aj_spl, inc_spl, Nsamples =get_aj_inc_star(stardir, confidence=confidence)
			aj_data[k][i,:,:]=aj_stats*unit_nHz
			inc_data[k][i,:]=inc_stats
			if i == 0:
				aj_samples.append(np.zeros((Ncombi[k], jmax, Nsamples)))
				inc_samples.append(np.zeros((Ncombi[k], Nsamples)))
				#a1ovGamma_samples.append(np.zeros((Ncombi[k], Nsamples)))
			aj_samples[k][i,:,:]=aj_spl*unit_nHz
			inc_samples[k][i,:]=inc_spl
			print('		(b) Extracting l=0 Frequencies...')
			freqs_med=get_freqs_median(stardir, use_synthese=True)
			nu_l0=freqs_med[:,0]
			print('		(c) Extracting samples of a1ovGamma...')
			loga1_ov_gamma_file=stardir+'/Widths/Samples_a1_ov_Widths.sav'
			a1ovGamma_s=get_a1ovGamma_atnumax(loga1_ov_gamma_file, nu_l0, numax_star, log_results=False, fast=True)
			#print('Fast mode: ', np.median(a1ovGamma_s))
			#a1ovGamma_s=get_a1ovGamma_atnumax(loga1_ov_gamma_file, nu_l0, numax_star, log_results=False, fast=False)
			#print('Slow mode: ', np.median(a1ovGamma_s))
			#a1ovGamma_samples[k][i,:]=a1ovGamma_s
			a1ovGamma_data[k][i]=np.median(a1ovGamma_s)*0.98
	
	print('C. Plots of bias...')
	inc_true=np.asarray(lookup_combi('i', combi_data[0], param_names[0]), dtype=float)
	posOK=np.where(inc_true >= inc0)
	Ninc=len(inc_true)
	yr=[0.1,0.8] # a1/Gamma yrange is for the moment fixed. But eventually, it should be determined by the a1ovGamma_true
	do_j=[0,1,3] # LIST OF THE aj+1 THAT WE PROCESS: j=0 is for a1, etc...
	coeff_circles=800#
	for j in do_j: # loop over the aj terms
		aj_true, Gamma_at_numax_true, a1ovGamma_true=get_true_vals(j+1, combi_data[0], param_names[0], unit_nHz) # This is to get aj_true. We do not care about 
		print(" j=", j+1, ":")		
		fig ,ax= plt.subplots(tight_layout=True)
		ax.set_xlabel('Inclination (deg)')
		ax.set_ylabel(r'$a_1/\Gamma'+ '$ (no unit)')
		#ax.set_ylabel('$a_'+str(j+1) + '$ (nHz)')
		ax.set_ylim(yr[0], yr[1])
		ax.set_xlim(-1, 91)
		ax.text(0.95, 0.92, r'$a_'+str(j+1) + '$= {0:0.2f} nHz'.format(aj_true[0]) , verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=14)#, fontweight='bold')
		minfill=min(yr)
		plt.fill_between([0, inc0], [yr[1],yr[1]], minfill, color='silver') # Add a gray area showing the discarded zone
		# Evaluate the color ranges for the plots 
		min_b=99999
		max_b=-99999
		min_b_true=99999
		max_b_true=-99999
		for k  in range(Nscenario):
			aj_true, Gamma_at_numax_true, a1ovGamma_true=get_true_vals(j+1, combi_data[k], param_names[k], unit_nHz)
			bias_true=aj_data[k][:,j,2]-aj_true[:]
			bias=bias_true # This is not necessarily the true bias depending on tricks used to enhance the visibility of the plot within the range of interest of inc
			if saturate_colors == True:
				print(" Saturation of colors requested. We will impose an upper limit on the bias based on values of bias above inc = inc0 =", inc0)
				#print("bias before saturation:", bias)
				bias=saturate_bias(bias, inc_true, inc0, extend_limit=0.) # Trick to focus the color range in the bias for the region of interest in terms of inclination (inc>inc0 <=> inc>30)
				#print("bias after saturation :", bias)
			if abs_err == False:
				if min_b > bias.min():
					min_b=bias.min()
				if min_b_true > bias_true.min():
					min_b_true=bias_true.min()
				if max_b < bias.max():
					max_b=bias.max()
				if max_b_true < bias_true.max():
					max_b_true=bias_true.max()	
			else:
				min_b=0
				min_b_true=0
				if max_b < np.abs(bias).max():
					max_b=np.abs(bias).max()
				if max_b_true < np.abs(bias_true).max():
					max_b_true=np.abs(bias_true).max()					
		#
		if abs_err == False:
			if min_b<0 and max_b>0:
				normalize = mcolors.TwoSlopeNorm(vcenter=0, vmin=min_b, vmax=max_b)
			if min_b>0 and max_b>0:
				normalize = mcolors.TwoSlopeNorm(vcenter=(min_b+max_b)/2, vmin=min_b, vmax=max_b)
			if min_b<0 and max_b<0:
				normalize = mcolors.TwoSlopeNorm(vcenter=(min_b+max_b)/2, vmin=min_b, vmax=max_b)
		else:
			normalize=mcolors.Normalize(vmin=0, vmax=max_b)
		#
		# Make the plots
		for k in range(Nscenario):
			aj_true, Gamma_at_numax_true, a1ovGamma_true=get_true_vals(j+1, combi_data[k], param_names[k], unit_nHz)
			#
			# Compute the mean error over inc>=inc0 and its standard deviation to get the mean expected uncertainty for each k
			zerr=make_error_from_stats(aj_data[k][:,j,:]) # The error on aj
			ax.text(79, a1ovGamma_true[0]-0.035, r"$\sigma = {0:0.0f} \pm {1:0.0f}$".format(np.median(zerr[:,np.where(inc_true>=inc0)]), np.std(zerr[:,np.where(inc_true>=inc0)])/np.sqrt(len(zerr[:,np.where(inc_true>=inc0)]))), fontsize=8)
			ax.plot(inc_true, a1ovGamma_true, linestyle='dashed',dashes=(5, 10), color='black') # Horizontal line
			bias=aj_data[k][:,j,2]-aj_true[:]
			if saturate_colors == True:
				bias=saturate_bias(bias, inc_true, inc0, extend_limit=0.) # Trick to focus the color range in the bias for the region of interest in terms of inclination (inc>inc0 <=> inc>30)
			if abs_err==False:
				cols=(bias - min_b)/(max_b - min_b)  # Normalise the range between 0 and 1
				colormap=ListedColormap(sns.diverging_palette(260, 0, s=150, l=50, n=5*len(cols), center='dark').as_hex())
				s=sns.scatterplot(x=inc_data[k][:,2], y=a1ovGamma_data[k][:], c=bias, norm=normalize, cmap=colormap, ax=ax, data=bias, s=coeff_circles*np.abs(bias/(max_b - min_b))) # Circles of size proportional to the bias
				scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
				scalarmappaple.set_array(bias)
			else:
				cols=(np.abs(bias) - min_b)/(max_b - min_b)  # Normalise the range between 0 and 1
				cols_true=(np.abs(bias_true) - min_b_true)/(max_b_true - min_b_true)  # Normalise the range between 0 and 1
				print("cols =", cols)
				print("cols_true =", cols_true)
				colormap=ListedColormap(sns.dark_palette(color='red').as_hex())
				#s=sns.scatterplot(x=inc_data[k][:,2], y=a1ovGamma_data[k][:], c=np.abs(bias), norm=normalize, cmap=colormap, ax=ax, data=np.abs(bias), s=coeff_circles*np.abs(bias_true/(max_b_true - min_b_true))) # Circles of size proportional to the bias
				s=sns.scatterplot(x=inc_data[k][:,2], y=a1ovGamma_data[k][:], c=np.abs(bias), norm=normalize, cmap=colormap, ax=ax, data=np.abs(bias), s=coeff_circles*cols) # Circles of size proportional to the bias
				# Specific color for the inclinations with known really bad estimates
				posBAD=np.where(inc_true<inc0)
				inc_bad=inc_true[posBAD] 
				for i in range(len(inc_bad)):
					#ax.scatter(inc_data[k][i,2], a1ovGamma_data[k][i], color=colorBAD, s=coeff_circles*np.abs(bias_true[i]/(max_b_true - min_b_true)))
					ax.scatter(inc_data[k][i,2], a1ovGamma_data[k][i], color=colorBAD, s=coeff_circles*cols[i])
				scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
				scalarmappaple.set_array(np.abs(bias))
			for i in range(Ninc):
				if inc_true[i]>=inc0:
					ax.annotate("", xy=(inc_data[k][i,2], a1ovGamma_data[k][i]), xytext=(inc_true[i], a1ovGamma_true[i]), color=colormap(cols[i]), va="center", ha="center", arrowprops=dict(arrowstyle="-", color=colormap(cols[i])))
				else:
					ax.annotate("", xy=(inc_data[k][i,2], a1ovGamma_data[k][i]), xytext=(inc_true[i], a1ovGamma_true[i]), color=colorBAD, va="center", ha="center", arrowprops=dict(arrowstyle="-", color=colorBAD))
				if bias_true[i] > 0: 
					ax.scatter(inc_data[k][i,2], a1ovGamma_data[k][i], marker='.', color='black')
				if bias_true[i] < 0:
					#ax.scatter(inc_data[k][i,2], a1ovGamma_data[k][i], marker='x', color='black', s=0.2*coeff_circles*np.abs(bias_true[i]/(max_b_true - min_b_true)))
					ax.scatter(inc_data[k][i,2], a1ovGamma_data[k][i], marker='x', color='black', s=0.2*coeff_circles*cols[i])
		fig.colorbar(scalarmappaple)
		#plt.show()
		plt.savefig(fileout+'_a'+str(j+1)+'.png', dpi=300)
		plt.close('all')
		print('File saved: ', fileout+'_a'+str(j+1)+'.png')

#dir_root='/Users/obenomar/tmp/test_a2AR/tmp/data/'
#combi_files=[dir_root +'HNR20_a1ovGamma0.5_Tobs730_Polar/Sources/Combinations.txt', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Equatorial/Sources/Combinations.txt']

## Looking at the Acquire phase
#MCMCdir=[dir_root+'/HNR20_a1ovGamma0.5_Tobs730_Polar/postMCMC/Level1/', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Equatorial/postMCMC/Level1/']
#fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/data/Result_Summary/All_plot'
##fileout='/Users/obenomar/tmp/test_a2AR/tmp/data/HNR20_a1ovGamma0.5_Tobs730_Equatorial/postMCMC/Level2/A/plot'

# --- Execution of the v1 -----
#j=0
#bias_analysis_v1(MCMCdir[j], combi_file[j], fileout=fileout)

# --- Execution of the v2 -----
#col_arrows=['cornflowerblue', 'deepskyblue']
#labels=['Polar cap', 'Equatorial Band']
#bias_analysis_v2(MCMCdir, combi_files, labels=labels, fileout=fileout_all, col_arrows=col_arrows)


# ----------------------------
# --- Execution of the v3 ----
# ----------------------------
dir_root='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/'
numax_star=2150.

## ----- TESTS ----
#combi_files=[dir_root +'HNR20_a1ovGamma0.4_Tobs730_Polar/Combinations.txt', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Polar/Combinations.txt', dir_root+'HNR20_a1ovGamma0.6_Tobs730_Polar/Combinations.txt']
## Looking at the Acquire phase
#MCMCdir=[dir_root+'/HNR20_a1ovGamma0.4_Tobs730_Polar/postMCMC/Level1/', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Polar/postMCMC/Level1/',dir_root+'HNR20_a1ovGamma0.6_Tobs730_Polar/postMCMC/Level1/']
#fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/data/Result_Summary/Bias_map'
#bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)
## --------

# -------- HNR 30 Tobs=730 days Polar -------
combi_files=[dir_root +'HNR30_a1ovGamma0.4_Tobs730_Polar/Combinations.txt', dir_root+'HNR30_a1ovGamma0.5_Tobs730_Polar/Combinations.txt', dir_root+'HNR30_a1ovGamma0.6_Tobs730_Polar/Combinations.txt']
MCMCdir=[dir_root+'/HNR30_a1ovGamma0.4_Tobs730_Polar/products/', dir_root+'HNR30_a1ovGamma0.5_Tobs730_Polar/products/',dir_root+'HNR30_a1ovGamma0.6_Tobs730_Polar/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR30_Tobs730_Polar'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)

# -------- HNR 30 Tobs=730 days Equatorial -------
combi_files=[dir_root +'HNR30_a1ovGamma0.4_Tobs730_Equatorial/Combinations.txt', dir_root+'HNR30_a1ovGamma0.5_Tobs730_Equatorial/Combinations.txt', dir_root+'HNR30_a1ovGamma0.6_Tobs730_Equatorial/Combinations.txt']
MCMCdir=[dir_root+'/HNR30_a1ovGamma0.4_Tobs730_Equatorial/products/', dir_root+'HNR30_a1ovGamma0.5_Tobs730_Equatorial/products/',dir_root+'HNR30_a1ovGamma0.6_Tobs730_Equatorial/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR30_Tobs730_Equatorial'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)

'''
# -------- HNR 20 Tobs=730 days Polar -------
combi_files=[dir_root +'HNR20_a1ovGamma0.4_Tobs730_Polar/Combinations.txt', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Polar/Combinations.txt', dir_root+'HNR20_a1ovGamma0.6_Tobs730_Polar/Combinations.txt']
MCMCdir=[dir_root+'/HNR20_a1ovGamma0.4_Tobs730_Polar/products/', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Polar/products/',dir_root+'HNR20_a1ovGamma0.6_Tobs730_Polar/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR20_Tobs730_Polar'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)

# -------- HNR 20 Tobs=730 days Equatorial -------
combi_files=[dir_root +'HNR20_a1ovGamma0.4_Tobs730_Equatorial/Combinations.txt', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Equatorial/Combinations.txt', dir_root+'HNR20_a1ovGamma0.6_Tobs730_Equatorial/Combinations.txt']
MCMCdir=[dir_root+'/HNR20_a1ovGamma0.4_Tobs730_Equatorial/products/', dir_root+'HNR20_a1ovGamma0.5_Tobs730_Equatorial/products/',dir_root+'HNR20_a1ovGamma0.6_Tobs730_Equatorial/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR20_Tobs730_Equatorial'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)

# -------- HNR 10 Tobs=730 days Polar -------
combi_files=[dir_root +'HNR10_a1ovGamma0.4_Tobs730_Polar/Combinations.txt', dir_root+'HNR10_a1ovGamma0.5_Tobs730_Polar/Combinations.txt', dir_root+'HNR10_a1ovGamma0.6_Tobs730_Polar/Combinations.txt']
MCMCdir=[dir_root+'/HNR10_a1ovGamma0.4_Tobs730_Polar/products/', dir_root+'HNR10_a1ovGamma0.5_Tobs730_Polar/products/',dir_root+'HNR10_a1ovGamma0.6_Tobs730_Polar/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR10_Tobs730_Polar'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)

# -------- HNR 10 Tobs=730 days Equatorial -------
combi_files=[dir_root +'HNR10_a1ovGamma0.4_Tobs730_Equatorial/Combinations.txt', dir_root+'HNR10_a1ovGamma0.5_Tobs730_Equatorial/Combinations.txt', dir_root+'HNR10_a1ovGamma0.6_Tobs730_Equatorial/Combinations.txt']
MCMCdir=[dir_root+'/HNR10_a1ovGamma0.4_Tobs730_Equatorial/products/', dir_root+'HNR10_a1ovGamma0.5_Tobs730_Equatorial/products/',dir_root+'HNR10_a1ovGamma0.6_Tobs730_Equatorial/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR10_Tobs730_Equatorial'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=None)
'''
# -------------
# -------------

# -------- HNR 30 Tobs=1460 days Polar -------
combi_files=[dir_root +'HNR30_a1ovGamma0.4_Tobs1460_Polar/Combinations.txt', dir_root+'HNR30_a1ovGamma0.5_Tobs1460_Polar/Combinations.txt', dir_root+'HNR30_a1ovGamma0.6_Tobs1460_Polar/Combinations.txt']
MCMCdir=[dir_root+'/HNR30_a1ovGamma0.4_Tobs1460_Polar/products/', dir_root+'HNR30_a1ovGamma0.5_Tobs1460_Polar/products/',dir_root+'HNR30_a1ovGamma0.6_Tobs1460_Polar/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR30_Tobs1460_Polar'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True)

# -------- HNR 30 Tobs=1460 days Equatorial -------
combi_files=[dir_root +'HNR30_a1ovGamma0.4_Tobs1460_Equatorial/Combinations.txt', dir_root+'HNR30_a1ovGamma0.5_Tobs1460_Equatorial/Combinations.txt', dir_root+'HNR30_a1ovGamma0.6_Tobs1460_Equatorial/Combinations.txt']
MCMCdir=[dir_root+'/HNR30_a1ovGamma0.4_Tobs1460_Equatorial/products/', dir_root+'HNR30_a1ovGamma0.5_Tobs1460_Equatorial/products/',dir_root+'HNR30_a1ovGamma0.6_Tobs1460_Equatorial/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR30_Tobs1460_Equatorial'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True)

'''
# -------- HNR 20 Tobs=1460 days Polar -------
combi_files=[dir_root +'HNR1020_a1ovGamma0.4_Tobs1460_Polar/Combinations.txt', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Polar/Combinations.txt', dir_root+'HNR20_a1ovGamma0.6_Tobs1460_Polar/Combinations.txt']
MCMCdir=[dir_root+'/HNR1020_a1ovGamma0.4_Tobs1460_Polar/products/', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Polar/products/',dir_root+'HNR20_a1ovGamma0.6_Tobs1460_Polar/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR20_Tobs1460_Polar'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=20)

# -------- HNR 20 Tobs=1460 days Equatorial -------
combi_files=[dir_root +'HNR1020_a1ovGamma0.4_Tobs1460_Equatorial/Combinations.txt', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Equatorial/Combinations.txt', dir_root+'HNR20_a1ovGamma0.6_Tobs1460_Equatorial/Combinations.txt']
MCMCdir=[dir_root+'/HNR1020_a1ovGamma0.4_Tobs1460_Equatorial/products/', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Equatorial/products/',dir_root+'HNR20_a1ovGamma0.6_Tobs1460_Equatorial/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR20_Tobs1460_Equatorial'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=20)

# -------- HNR 10 Tobs=1460 days Polar -------
combi_files=[dir_root +'HNR1020_a1ovGamma0.4_Tobs1460_Polar/Combinations.txt', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Polar/Combinations.txt', dir_root+'HNR10_a1ovGamma0.6_Tobs1460_Polar/Combinations.txt']
MCMCdir=[dir_root+'/HNR1020_a1ovGamma0.4_Tobs1460_Polar/products/', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Polar/products/',dir_root+'HNR10_a1ovGamma0.6_Tobs1460_Polar/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR10_Tobs1460_Polar'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=10)

# -------- HNR 10 Tobs=1460 days Equatorial -------
combi_files=[dir_root +'HNR1020_a1ovGamma0.4_Tobs1460_Equatorial/Combinations.txt', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Equatorial/Combinations.txt', dir_root+'HNR10_a1ovGamma0.6_Tobs1460_Equatorial/Combinations.txt']
MCMCdir=[dir_root+'/HNR1020_a1ovGamma0.4_Tobs1460_Equatorial/products/', dir_root+'HNR1020_a1ovGamma0.5_Tobs1460_Equatorial/products/',dir_root+'HNR20_a1ovGamma0.6_Tobs1460_Equatorial/products/']
fileout_all='/Users/obenomar/tmp/test_a2AR/tmp/Simulationdata/Result_Summary/Bias_map_HNR10_Tobs1460_Equatorial'
bias_analysis_v3(MCMCdir, combi_files, numax_star, labels=None, fileout=fileout_all, abs_err=True, saturate_colors=True, filter_HNR=10)
'''
exit()
