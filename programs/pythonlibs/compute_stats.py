# General set of functions generating statistical summaries
import numpy as np

def make_stats(samples, confidence=[2.25,16,50,84,97.75]):
	N=len(samples)
	s=np.sort(samples)
	cdf = 100.*np.array(range(N))/float(N) # in %
	r=np.interp(confidence, cdf, s)
	return r

def make_error_from_stats(stats):
	err=np.zeros((2, len(stats[:,0]))) 
	err[0,:]=stats[:, 2] - stats[:, 1]
	err[1,:]=stats[:, 3] - stats[:, 2]
	return err

def uncertainties(data):
    try:
        Nparams=len(data[0,:])
        errors=np.zeros((2,Nparams))
        med=np.zeros(Nparams)
        stats_all=np.zeros((Nparams,5))
        for i in range(Nparams):
            print(data[:,i])
            stats = np.percentile(data[:, i], [2.25, 16, 50, 84, 97.75])
            errors[0,i]=stats[2] - stats[1]
            errors[1,i]=stats[3] - stats[2]
            med[i]=stats[2]
            stats_all[i,:]=stats  
            return med, errors, stats_all  
    except:
        errors=np.zeros(2)
        stats = np.percentile(data, [2.25, 16, 50, 84, 97.75])
        errors[0]=stats[2] - stats[1]
        errors[1]=stats[3] - stats[2]
        return stats[2], errors, stats 