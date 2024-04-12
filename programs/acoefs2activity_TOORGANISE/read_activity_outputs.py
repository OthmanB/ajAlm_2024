'''
NOTE: DO NOT CALL THIS FUNCTION. CALL INSTEAD read_activity_outputs_emcee.py
      THIS FUNCTION IS INTENDED TO BE ERASED AT SOME POINT
'''

import numpy as np

def read_logposterior(file_name='logposterior.npy'):
	'''
		Read the logposterior file and return all of the logposterior probabilities
	'''
	lpdf=np.load(file_name)
	return lpdf


def read_params(file, relax=np.asarray([True, True, False, True])):
	'''
		Use read_mcmc_posteriors and shape them in a suitable way for MLCJ.py
	'''
	smcmc, labels, Nparams, Nsamples=read_mcmc_posteriors(file)
	sinput=smcmc[:,0] # recover the initial values... this will be used to construct pref_all
	pvars=np.where(relax == True)
	smcmc=smcmc[np.reshape(pvars, np.size(pvars)), :] # keep only the variables (save memory...)
	return smcmc, sinput, pvars

def read_mcmc_posteriors(file):
	'''
		Reads the samples of the posteriors for each parameters
	'''
	smcmc=np.load(file)
	Nparams=len(smcmc[0,:])
	Nsamples=len(smcmc[:,0])
	labels = ["epsilon_nl0", "epsilon_nl1", "theta0", "delta"]
	return smcmc, labels, Nparams, Nsamples	


def get_ref_MAP(lpdf, smcmc):
	'''
		Returns the maximum a posteriori and the associated values of the parameters, provided:
			- lpdf: The logposterior probability of all of the samples
			- smcmc: The array of parameters used during the fit
	'''
	pos= np.unravel_index(lpdf.argmax(), lpdf.shape) # get the position of the maximum of the lpdf
	#print("pos=", pos)
	#print("shape(smcmc): ", np.shape(smcmc))
	return lpdf[pos[0]], smcmc[:,pos[0]], pos

def get_ref_nearmean_lpdf(lpdf, smcmc):
	'''
		Return the value of the mean of the logposterior and the nearest value of the parameters for that mean
			- lpdf: The logposterior probability of all of the samples
			- smcmc: The array of parameters used during the fit
	'''
	lpdf_mean=np.mean(lpdf)
	delta=np.abs(lpdf-lpdf_mean)
	pos= np.unravel_index(delta.argmin(), delta.shape) # get the position of the minimum of the difference betweeb lpdf_mean and lpdf
	return lpdf[pos[0]], smcmc[:,pos[0]], pos[0]

def compute_covarmat(smcmc):
	'''
		Use the samples to estimate the covariance matrix for the variables
	'''
	c=np.cov(smcmc) # VERIFY THAT THERE IS NO TRANSPOSE ISSUE
	return c

def reduce_samples(smcmc, lpdf):
# look at the occurence of lpdf values (due to the fact that ~25% only are accepted)
# and return the occurence rate of each values along with reduced vectors smcmc_reduced and lpdf_reduced
	lpdf_tmp=lpdf
	smcmc_tmp=smcmc
	#
	lpdf_reduced=np.empty(len(lpdf))
	smcmc_reduced=np.empty([len(smcmc[:,0]),len(smcmc[0,:])])
	recur=np.empty(len(lpdf))
	#
	i=0
	print("reducing arrays...")
	while len(lpdf_tmp) != 0: # process the elements
		pos_r=np.where(lpdf_tmp == lpdf_tmp[0])
		r=len(pos_r[0])
		lpdf_reduced[i]=lpdf_tmp[0]
		smcmc_reduced[i,:]=smcmc_tmp[0,:]
		recur[i]=r
		#
		lpdf_tmp=np.delete(lpdf_tmp, pos_r)
		smcmc_tmp=np.delete(smcmc_tmp, pos_r, axis=0)
		#
		#new_pos=np.where(lpdf_tmp != lpdf_tmp[0])
		#lpdf_tmp=lpdf_tmp[new_pos]
		#smcmc_tmp=smcmc_tmp[new_pos,:]
		#print "pos_r=", pos_r[0], "  r=", r
		if i%1000 == 0:
			print("i=",i)
		i=i+1
	#
	pos_del=np.arange(len(lpdf) - i) + i
	lpdf_reduced=np.delete(lpdf_reduced, pos_del)
	smcmc_reduced=np.delete(smcmc_reduced, pos_del, axis=0)
	recur=np.delete(recur, pos_del)
	#print("last element of lpdf_reduced ", lpdf_reduced[i-1])
	#
	return np.array(smcmc_reduced), np.array(lpdf_reduced), np.array(recur)

