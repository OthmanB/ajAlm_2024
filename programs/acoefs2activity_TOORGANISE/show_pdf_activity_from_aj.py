'''
	This is the ensemble of function that are used to
	read MCMC files from the EMCMC run (fit_ajsig.py program) 
	and to make paper-ready plots for the activity parameters
	epsilon, theta0 and delta
	At the end of this file, there is a list of calls for the case of:
		- The Sun at minimum and maximum of activity
		- 16 Cyg B for the two solutions of a4
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import interpolate
from scipy.signal import medfilt
from  matplotlib import ticker
from nice_hist import nice_hist
from read_activity_outputs_emcee import read_mcmc_posteriors
from read_outputs_tamcmc import read_params_tamcmc
from scipy.signal import savgol_filter
from fit_a2sig import prior_uniform
from fit_a2sig import prior_jeffreys

def filter_posteriors(data, labels, Nparams, Nsamples, burnin, skip_epsilon_nl1=False):#, keep_positive=False, keep_epsilon=None):
	Nburnin=int(burnin*Nsamples) # USER CAN REMOVE A FRACTION OF THE SAMPLES AT THE BEGINING
	if burnin > 1:
		print("Error in filtering posteriors: You specified a value of the burnin parameteter exceeding 1")
		print("                               This should be a fraction of Nsamples and must be < 1")
		exit()
	# -- Section to skip epsilon_nl1 --
	if skip_epsilon_nl1 == True:
		data1=np.zeros((Nsamples-Nburnin, Nparams-1))
		data1[:,0]=data[Nburnin:,0]
		data1[:,1]=data[Nburnin:,2]
		data1[:,2]=data[Nburnin:,3]
		Nparams=len(data1[0,:])	
	else:
		data1=np.zeros((Nsamples-Nburnin, Nparams))
		for i in range(Nparams):
			data1[:,i]=data[Nburnin:,i]
	labels=["epsilon", "theta", "delta"]
	return data1, labels, Nparams, Nsamples

def get_uncertainties(data, Nparams, labels, rad2deg=True):
	errors=np.zeros((2,Nparams))
	med=np.zeros(Nparams)
	stats_all=np.zeros((Nparams,5))
	for i in range(Nparams):
		stats = np.percentile(data[:, i], [2.25, 16, 50, 84, 97.75])
		print(' stats[',i,'] = ', stats)
		if labels[i] == 'epsilon':
			errors[0,i]=stats[2] - stats[1]
			errors[1,i]=stats[3] - stats[2]
			med[i]=stats[2]
			stats_all[i,:]=stats
		else:
			if rad2deg == True:
				errors[0,i]=(stats[2] - stats[1])*180./np.pi
				errors[1,i]=(stats[3] - stats[2])*180./np.pi
				med[i]=stats[2]*180./np.pi
				stats_all[i,:]=stats*180./np.pi
			else:
				errors[0,i]=stats[2] - stats[1]
				errors[1,i]=stats[3] - stats[2]
				med[i]=stats[2]
				stats_all[i,:]=stats
	return med, errors, stats_all


def show_pdf_vpaper(dir_mcmc=['/Users/obenomar/tmp/test_a2AR/tmp/results_activity/gate_8379927_a2a4_fast/'], names=[None], file_out='plots', burnin=0, do_inset=True, logepsilon=[False, False], xlim_epsilon=[[0, 1e-2],[0, 1e-2]], text_index=None, data_format=['emcee']):
	'''
		Main function to make the pdfs plots for the paper Benomar2022+ and for the activity parameters: epsilon_nl, theta0, delta
		dir_mcmc: Array with directories that specify where to look at the data. 
		file_out: Name of the output image file
		burnin: Fraction (between 0 and 1) of data points that are discarded at the begning of the serie of samples. Allows to remove some burnin, if necessary
		logepsilon[0]: Show the main plot of epsilon in log, if set to True
		logepsilon[1]: Show the main plot of epsilon in log, if set to True
	'''
	if names[0] == None:
		names=np.linspace(1, len(dir_mcmc), len(dir_mcmc))
	if len(data_format) > 7: # We have the prior pdfs as well in that case
		data_priors, sinput,labels=read_params_tamcmc(data_format[6], data_format[7], phase=data_format[8], chain=data_format[9], first_index=data_format[10], period=data_format[11], epsi_iscte=False)
	#
	data_s=[]
	labels_s=[]
	Nparams_s=[]
	Nsamples_s=[]
	med_s=[]
	errors=[]
	statistics=[]
	for j in range(len(dir_mcmc)):
		d=dir_mcmc[j]
		print("    - - - - - - dir_mcmc : - - - - - - - ")
		print("        d=", d)
		if data_format[0] == 'emcee':
			data, labels, Nparams, Nsamples=read_mcmc_posteriors(d + 'samples.npy')
			rad2deg=True # The emcee processes are giving theta and delta in rad
		else:
			# If 'tamcmc', the expected format is: process_name= data_format[1], phase= data_format[2], chain = data_format[3], etc...
			data, sinput,labels=read_params_tamcmc(d, data_format[1][j], phase=data_format[2], chain=data_format[3], first_index=data_format[4], period=data_format[5], epsi_iscte=False)
			Nparams=len(data[0,:])
			Nsamples=len(data[:,0])
			rad2deg=False # The TAMCMC processes are already giving theta and delta in degrees
			#print(Nsamples)
			#exit()
		if Nparams == 4:
			data, labels, Nparams,Nsamples=filter_posteriors(data, labels, Nparams,Nsamples, burnin, skip_epsilon_nl1=True)
		else:
			data, labels, Nparams,Nsamples=filter_posteriors(data, labels, Nparams,Nsamples, burnin, skip_epsilon_nl1=False)
		#
		# Evaluate uncertainties using the samples
		med,err, stats=get_uncertainties(data, Nparams, labels, rad2deg=rad2deg)
		#
		# Show summary on the statistics
		string='#param   median   err_m     err_p\n'
		for i in range(Nparams):
			print(labels[i], ' =', med[i], '  (-  ', err[0, i], '  ,  +', err[1,i], ')')
			string=string + '  {}  {:0.6f}      {:0.6f}      {:3.6f}\n'.format(labels[i], med[i], err[0,i], err[1,i])
		fsum=open(file_out+'_'+names[j]+'_summary.txt', 'w')
		fsum.write(string)
		fsum.close()
		#
		data_s.append(data)
		labels_s.append(labels)
		Nparams_s.append(Nparams)
		Nsamples_s.append(Nsamples)
		med_s.append(med)
		errors.append(err)
		statistics.append(stats)
	#
	cols=[['black', 'gray', 'darkgray'], ['red', 'tomato','lightsalmon']]
	alphas=[None, 0.4]
	bins=20#75
	for i in range(Nparams):
		fig_1d, ax = plt.subplots(1, figsize=(12, 6))
		fig_bin, ax_bin = plt.subplots(1, figsize=(12, 6))
		if len(data_format) > 7: # Handle the case with prior given
			h, x=np.histogram(data_s[0][:,i], bins=70*bins, density=True)
			xbins=np.linspace(0, np.max(data_priors[:,i]),70*bins)
			hist_priors, bins_prior=np.histogram(data_priors[:,i], bins=xbins, density=True)
			scale=max(h)/max(hist_priors)
			shist = savgol_filter(hist_priors, 40*bins+1, 3)
			ax.plot(xbins[0:-1], shist*scale, linestyle='--', color='gray')
			#y,x,p=ax.hist(data_priors[:,i], linestyle='--', bins=70*bins, histtype='step', color='gray', density=True)
		for j in range(len(dir_mcmc)):
			if text_index != None:
				ax.annotate(text_index[i], xy=(0.05, 0.92), xycoords=ax.transAxes, fontsize=22)
			if labels_s[j][i] == "epsilon":
				ax.set_xlim([xlim_epsilon[0][0], xlim_epsilon[0][1]])
				if j == 0 and do_inset == True:
					#ax_inset = inset_axes(ax, width="30%", height="40%", loc=1)
					ax_inset = inset_axes(ax, width="30%", height="40%", bbox_to_anchor=(-0.05,-0.05,1,1), bbox_transform=ax.transAxes) 
					ax_inset.set_yticks([])
					ax_inset.axhline(0, linestyle='--', color='black')
					ax_inset.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
				nice_hist(data_s[j][:,i], statistics[j][i,:], ax=ax, intervals=[True,False], binning=10*bins, color=cols[j],alpha=alphas[j])
				if do_inset == True:
					y,x,p=ax_bin.hist(data_s[j][:,i], linestyle='-', bins=70*bins, histtype='step', color=cols[j][0], density=True)
					ax_inset.step(x[0:len(y)],medfilt(y, 1), color=cols[j][0], linestyle='-')
					ax_inset.set_xlim([xlim_epsilon[1][0], xlim_epsilon[1][1]])
					#ax_inset.axvline(med_s[j][i], linestyle='--', color=cols[j][0])
					if logepsilon[1] == True:
						ax_inset.set_xscale('log')
						ax_inset.set_xlim([xlim_epsilon[1][0], xlim_epsilon[1][1]])
				if logepsilon[0] == True:
		 			ax.set_xscale('log')
			else:
				if labels_s[j][i] == "delta":
					Nbins=8*bins
				else:
					Nbins=2*bins
				if rad2deg == True:
					nice_hist(data_s[j][:,i]*180/np.pi, statistics[j][i,:], ax=ax, intervals=[True,False], binning=Nbins,color=cols[j], alpha=alphas[j])
				else:
					nice_hist(data_s[j][:,i], statistics[j][i,:], ax=ax, intervals=[True,False], binning=Nbins,color=cols[j], alpha=alphas[j])					
			if labels_s[j][i] == 'epsilon':
				ax.set_xlabel(r'$\epsilon$ (no unit)', fontsize=18)
			if labels_s[j][i] == 'theta':
				ax.set_xlabel(r'$\theta$ (deg)', fontsize=18)
			if labels_s[j][i] == 'delta':
				ax.set_xlabel(r'$\delta$ (deg)', fontsize=18)
			ax.set_ylabel("PDF", fontsize=20)
			ax.tick_params(axis='x', labelsize=18)
			ax.tick_params(axis='y', labelsize=18)		
			fig_1d.savefig(file_out+ '_' + names[0] + '_pdf_'+ str(i) + '.jpg')
	print('Saved files with syntax: ', file_out)
	

# ----------------------------------------
# --------------- USE OF TAMCMC -----------
# ----------------------------------------
# Common out directory
dir_out='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/PostMCMC/'

# The Sun 1999-2002
dir_mcmc=['/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/'] # Must be an array
#dir_mcmc_CFonly='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/' # Must be a single string
process_names=['aj_ARonly_Sun_19992002_l12only']
#process_name_CFonly='aj_CFonly_Sun_19992002_l12only' # Must be a single string
names=['19992002_incfix_fast_Priorevalrange_a2ARonly']
data_format=['tamcmc', process_names, 'A', 0, 0, 1]#, dir_mcmc_CFonly, process_name_CFonly, 'A', 0, 0, 1] # data_type, process_name, phase, chain, period + Priors info: dir_mcmc_CFonly, process_name, 'A', 0,0,1
file_out=dir_out+'plots'
show_pdf_vpaper(dir_mcmc=dir_mcmc,names=names, logepsilon=[False,False], burnin=0.0, xlim_epsilon=[[0, 5e-3],[0,5e-4]], text_index=['(a)', '(b)', '(c)'], data_format=data_format, file_out=file_out)

dir_mcmc=['/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/']
#dir_mcmc_CFonly='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/'
process_names=['aj_Sun-20062009_l12only']
#process_name_CFonly='aj_Sun-20062009_CFonly_l12only' # Must be a single string
names=['20062009_incfix_fast_a2ARonly']
data_format=['tamcmc', process_names, 'A', 0, 0, 1]#, dir_mcmc_CFonly, process_name_CFonly, 'A', 0, 0, 1] # data_type, process_name, phase, chain, period + Priors info: dir_mcmc_CFonly, process_name, 'A', 0,0,1
file_out=dir_out+'plots'
show_pdf_vpaper(dir_mcmc=dir_mcmc,names=names, logepsilon=[False,False], burnin=0.0, xlim_epsilon=[[0, 5e-3],[0,5e-4]], text_index=['(a)', '(b)', '(c)'], data_format=data_format, file_out=file_out)

#    16 Cyg A
dir_mcmc=['/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/', '/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/']
#dir_mcmc_CFonly='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/'
process_names=['aj_ARonly_16CygA', 'aj_ARonly_16CygA_bias']
#process_name_CFonly='aj_Sun-20062009_CFonly_l12only' # Must be a single string
names=['aj_ARonly_16CygA','aj_ARonly_16CygA_bias']
data_format=['tamcmc', process_names, 'A', 0, 0, 1]#, dir_mcmc_CFonly, process_name_CFonly, 'A', 0, 0, 1] # data_type, process_name, phase, chain, period + Priors info: dir_mcmc_CFonly, process_name, 'A', 0,0,1
file_out=dir_out+'plots'
show_pdf_vpaper(dir_mcmc=dir_mcmc,names=names, logepsilon=[False,False], burnin=0.0, xlim_epsilon=[[0, 5e-3],[0,5e-4]], text_index=['(a)', '(b)', '(c)'], data_format=data_format, file_out=file_out)
'''

#  16 Cyg B  LowSol
dir_mcmc=['/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/', '/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/']
#dir_mcmc_CFonly='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/'
process_names=['aj_ARonly_16CygB_LowSol', 'aj_ARonly_16CygB_LowSol_bias']
#process_name_CFonly='aj_Sun-20062009_CFonly_l12only' # Must be a single string
names=['aj_ARonly_16CygB_LowSol','aj_ARonly_16CygB_LowSol_bias']
data_format=['tamcmc', process_names, 'A', 0, 0, 1]#, dir_mcmc_CFonly, process_name_CFonly, 'A', 0, 0, 1] # data_type, process_name, phase, chain, period + Priors info: dir_mcmc_CFonly, process_name, 'A', 0,0,1
file_out=dir_out+'plots'
show_pdf_vpaper(dir_mcmc=dir_mcmc,names=names, logepsilon=[False,False], burnin=0.0, xlim_epsilon=[[0, 5e-3],[0,5e-4]], text_index=['(a)', '(b)', '(c)'], data_format=data_format, file_out=file_out)

#  16 Cyg B  UpSol
dir_mcmc=['/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/', '/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/']
#dir_mcmc_CFonly='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/TAMCMC/'
process_names=['aj_ARonly_16CygB_UpSol', 'aj_ARonly_16CygB_UpSol_bias']
#process_name_CFonly='aj_Sun-20062009_CFonly_l12only' # Must be a single string
names=['aj_ARonly_16CygB_UpSol','aj_ARonly_16CygB_UpSol_bias']
data_format=['tamcmc', process_names, 'A', 0, 0, 1]#, dir_mcmc_CFonly, process_name_CFonly, 'A', 0, 0, 1] # data_type, process_name, phase, chain, period + Priors info: dir_mcmc_CFonly, process_name, 'A', 0,0,1
file_out=dir_out+'plots'
show_pdf_vpaper(dir_mcmc=dir_mcmc,names=names, logepsilon=[False,False], burnin=0.0, xlim_epsilon=[[0, 5e-3],[0,5e-4]], text_index=['(a)', '(b)', '(c)'], data_format=data_format, file_out=file_out)
'''

exit()