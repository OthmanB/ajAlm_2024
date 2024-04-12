import numpy as np
from subprocess import Popen, PIPE
import os
import time

'''
Functions to directly read the outputs from a mcmc run
'''

def read_covarmat(dir_tamcmc_outputs, process_name, phase='A', chain=0):

	file2=dir_tamcmc_outputs + '/' + process_name + '/restore/' + process_name + '_restore_' + phase + '_2.dat'  # FILE WITH SIGMA
	file3=dir_tamcmc_outputs + '/' + process_name + '/restore/' + process_name + '_restore_' + phase + '_3.dat'  # File with the covarmat
	Nchains, Nvars, iteration, variable_names, sigmas, sigmas_mean, mus, mus_mean=read_restore_file2(file2)
	Nchains, Nvars, iteration, variable_names, covarmats=read_restore_file3(file3)
	c=sigmas[chain] * covarmats[chain, :, :] # Calculates the covariance matrix of the coldest chain using sigma and the covmat
	#c=covarmats[chain, :, :]
	return c

def read_restore_file2(filename):
	'''
		Read the restore files with the index = 2
		This file contains sigma and mu, the scaling factors of the MALA algorithm
	'''
	f=open(filename)
	txt=f.read()
	f.close()
	data=txt.split('\n')
	header=[]
	loglikelihood=[]
	logprior=[]
	logposterior=[]
	for d in data:
		s=d.split()
		if len(s) !=0:
			if s[0] == '#':
				header.append(d)
			if s[0] == '!':
				v=s[1].split('=') # Separate the value from the rest to be able to pick the value
				if v[0] == 'Nchains':
					Nchains=int(d.split('=')[-1]) 
				if v[0] == 'Nvars':
					Nvars=int(d.split('=')[-1]) 
				if v[0] == 'iteration':
					iteration=int(d.split('=')[-1]) 
				if v[0] == 'variable_names':
					st=d.split('=')[-1]
					variable_names=st.split()					
				if v[0] == 'sigmas':
					st=d.split('=')[-1]
					sigmas=np.asarray(st.split(), dtype=float)
				if v[0] == 'sigmas_mean':
					st=d.split('=')[-1]
					sigmas_mean=np.asarray(st.split(), dtype=float)
	for i in range(len(data)):
		s=data[i].split()
		if len(s) !=0:
			if s[0] == '!':
				v=s[1].split('=') # Separate the value from the rest to be able to pick the value
				if v[0] == 'mus':
					mus_raw=data[i+1:i+1+Nchains]
					mus=np.zeros((Nchains, Nvars))
					for k in range(len(mus_raw)):
						mus[k,:]=mus_raw[k].split()
				if v[0] == 'mus_mean':
					mus_raw=data[i+1:i+1+Nchains]
					mus_mean=np.zeros((Nchains, Nvars))
					for k in range(len(mus_raw)):
						mus_mean[k,:]=mus_raw[k].split()
	return Nchains, Nvars, iteration, variable_names, sigmas, sigmas_mean, mus, mus_mean

def read_restore_file3(filename):
	'''
		Read the restore files with the index = 3
		This file contains the covarmats of the MALA algorithm
	'''
	f=open(filename)
	txt=f.read()
	f.close()
	data=txt.split('\n')
	header=[]
	loglikelihood=[]
	logprior=[]
	logposterior=[]
	for d in data:
		s=d.split()
		if len(s) !=0:
			if s[0] == '#':
				header.append(d)
			if s[0] == '!':
				v=s[1].split('=') # Separate the value from the rest to be able to pick the value
				if v[0] == 'Nchains':
					Nchains=int(d.split('=')[-1]) 
				if v[0] == 'Nvars':
					Nvars=int(d.split('=')[-1]) 
				if v[0] == 'iteration':
					iteration=int(d.split('=')[-1]) 
				if v[0] == 'variable_names':
					st=d.split('=')[-1]
					variable_names=st.split()				
	for i in range(len(data)):
		s=data[i].split()
		if len(s) !=0:
			if s[0] == '!':
				v=s[1].split('=') # Separate the value from the rest to be able to pick the value
				if v[0] == 'covarmats':
					covarmats=np.zeros((Nchains, Nvars, Nvars))
					i0=i+2					
					for j in range(Nchains):
						covarmat_raw=data[i0:i0+Nvars]
						for k in range(Nvars):
							covarmats[j, k,:]=covarmat_raw[k].split()
						i0=i0+Nvars+1
	return Nchains, Nvars, iteration, variable_names, covarmats


def getstats_bin(dir_tamcmc_outputs, process_name, phase='A', chain=0, erase_tmp=True):
	'''
		Convert a binary output file from the tamcmc process into a simple text format
		The function makes use of the get_stats program created upon compilation of the tamcmc program
		dir_tamcmc_outputs: The roott directory that contains all of the outputs. The process outputs must be in a subdirectory of this directory
		core_filename: The name of the process that was ran (eg. 'aj_ARonly_Sun'). This is the name listed in the config_presets.cfg of the tamcmc program
		phase: The phase of the analysis
		chain: The chain index that need to be unpacked
		outfile: The name of the output file
	'''
	outfile='tmp/posteriors.txt'
	core_filename=dir_tamcmc_outputs + '/' + process_name + '/outputs/' + process_name + '_' + phase + '_stat_criteria' 
	process = Popen(["cpp_prg/./getstats", core_filename, str(chain), outfile], stdout=PIPE, stderr=PIPE)
	time.sleep(1)
	loglikelihood, logprior, logposterior, header, labels=getstats_txt(outfile)
	time.sleep(1)
	if erase_tmp == True:
		process = os.remove(outfile)
	print(header)
	print('--')
	return loglikelihood, logprior, logposterior

def bin2txt(dir_tamcmc_outputs, process_name, phase='A', chain=0, first_index=0, last_index=-1, period=1, erase_tmp=True, version='2022'):
	'''
		A function that calls bin2txt and read the outputs in order to return the samples for each parameters in a matrix form
	'''
	outdir='tmp/'
	core_filename=dir_tamcmc_outputs + '/' + process_name + '/outputs/' + process_name + '_' + phase + '_params'
	if version != "2022" : # The old version below 1.83.2 did not have 'last_index' as a parameter
		process = Popen(["cpp_prg/./bin2txt", core_filename, str(chain), outdir, str(first_index), str(period)], stdout=PIPE, stderr=PIPE)	
	else: 
		process = Popen(["cpp_prg/./bin2txt", core_filename, str(chain), outdir, str(first_index), str(last_index), str(period)], stdout=PIPE, stderr=PIPE)	
	
	time.sleep(1.)
	#
	# Iterate over all of the available files in the directory
	files=get_files_list(outdir, extension='.ASCII')
	if len(files) != 0:
		print(' bin2txt executed successfully')
	else:
		print('Error: No ASCII files in the tmp directory.')
		print('       It is likely that bin2txt failed to read the TAMCMC data')
		print('       Check manually the settings of bin2txt:')
		print('           - Program Location: cpp_prg/./bin2txt' )
		print('           - Root File Names : ', core_filename )
		print('           - Chain Index     : ', chain)
		print('           - Output Directory: ', outdir)
		print('           - First Index     : ', first_index)
		print('           - Period          : ', period)
		print('       Executed command: ')
		print('           cpp_prg/./bin2txt ', core_filename, ' ', str(chain), ' ', outdir, ' ', str(first_index), ' ', str(period))
		print(' --- ')	
		print('       Returned at Execution: ')
		res = process.communicate()
		message=str(res[0], "utf-8").split('\n')
		for m in message:
				print(m)
		exit()

	# Get the number of samples from the first file
	samples, varname=read_bin2txt_out(outdir + '000.ASCII')
	Nsamples=len(samples)
	if Nsamples == 1: # If the first parameter was fixed, try with the second one
		samples, varname=read_bin2txt_out(outdir + '001.ASCII')
		Nsamples=len(samples)
	print('Nsamples=', Nsamples)
	#
	smcmc=np.zeros((Nsamples, len(files)))
	labels=[]
	for f in files:
		index=int(f.split('.')[0]) # index of the parameter
		samples, varname=read_bin2txt_out(outdir + f)
		labels.append(varname)
		smcmc[:,index]=samples
	time.sleep(1)
	# erase files in the temporary directory
	if erase_tmp == True:
		for f in files:
			process = os.remove(outdir+f)
		process = os.remove(outdir+'plength.txt')
	return smcmc, labels

def test_bin2txt():
	dir_tamcmc_outputs='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/activity_TAMCMC/28-06-2022/ajfits/'
	process_name='aj_ARonly_Sun_19992002'
	bin2txt(dir_tamcmc_outputs, process_name, phase='A', chain=0, first_index=0, period=1, erase_tmp=True)

def read_bin2txt_out(file):
	# Function that get the header and the data within the text output files generated by bin2txt.cpp
	f=open(file)
	txt=f.read()
	f.close()
	data=txt.split('\n')
	varname=data[0].split('=')[1:]
	samples=[]
	for d in data[1:]:
		if d != '':
			samples.append(d)	
	return np.asarray(samples), varname

def getstats_txt(filename):
	'''
		A function that is able to read the text file created when running getstats (the tamcmc program that unpack data from the bin files)
	'''
	f=open(filename)
	txt=f.read()
	f.close()
	data=txt.split('\n')
	header=[]
	loglikelihood=[]
	logprior=[]
	logposterior=[]
	for d in data:
		s=d.split()
		if len(s) !=0:
			if s[0] == '#':
				header.append(d)
			else:
				loglikelihood.append(float(s[0]))
				logprior.append(float(s[1]))
				logposterior.append(float(s[2]))
	labels=header[-1]
	header=header[0:-2]
	loglikelihood=np.asarray(loglikelihood)
	logprior=np.asarray(logprior)
	logposterior=np.asarray(logposterior)
	return loglikelihood, logprior, logposterior, header, labels


def get_files_list(rootdir, extension='.ASCII'):
	'''
		A tiny function that scan a directory in order to find all of 
		its files
	'''
	files=[]
	for x in os.listdir(rootdir):
		if x.endswith(extension):
			files.append(x)
	return files


##########

def read_logposterior_tamcmc(dir_tamcmc_outputs, process_name, phase='A', chain=0):
	'''
		Read the logposterior file and return all of the logposterior probabilities
		Uses the getstats_bin() function to extract the data, convert them into text,
		read them with getstats_txt() and recover only what is relevant for MLCJ, ie
		the logposteriors (lpdf)
		dir_tamcmc_outputs: The roott directory that contains all of the outputs. The process outputs must be in a subdirectory of this directory
		core_filename: The name of the process that was ran (eg. 'aj_ARonly_Sun'). This is the name listed in the config_presets.cfg of the tamcmc program
		phase: The phase of the analysis
		chain: The chain index that need to be unpacked
		outfile: The name of the output file
	'''
	ll, lp, lpdf=getstats_bin(dir_tamcmc_outputs, process_name, phase=phase, chain=chain)
	return lpdf

def read_params_tamcmc(dir_tamcmc_outputs, process_name, phase='A', chain=0, first_index=0, period=1, epsi_iscte=False):
	'''
		Use bin2txt and shape them in a suitable way for MLCJ.py
		This function basically removes all of the constant terms.
	'''
	s, l=bin2txt(dir_tamcmc_outputs, process_name, phase=phase, chain=chain, first_index=first_index, period=period, erase_tmp=True)
	# --- Keeps only epsilon_0 (if epsi_iscte == False), theta0 and delta ---
	if epsi_iscte == False:
		smcmc=s[:,0:3] # keep epsilon_0, theta0, delta
		labels=l[0:3]
	else:
		smcmc=s[:,1:3] # keep theta0 and delta
		labels=l[1:3]
	sinput=smcmc[0,:] # recover the initial values... this will be used to construct pref_all
	return smcmc, sinput,labels


