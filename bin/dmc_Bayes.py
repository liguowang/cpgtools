#!/usr/bin/env python
'''
#=========================================================================================
Description
-----------
Different from statistical testing, this program tries to estimates "how different the
means between the two groups are" using Bayesian approach. An MCMC is used to estimate the
"means", "difference of means", "95% HDI (highest posterior density interval)", and the
osterior probability that the HDI does NOT include "0". 

It is similar to John Kruschke's BEST algorithm (Bayesian Estimation Supersedes T test)
(http://www.indiana.edu/~kruschke/BEST/).

Notes
-----
This program is much slower than T test due to MCMC (Markov chain Monte Carlo) step.
Running it with multiple threads is highly recommended.

'''
import sys,os
import subprocess
import numpy as np
from scipy import stats
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED
from cpgmodule import padjust
from multiprocessing import Process, Manager, current_process

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def dt(x, mu, sig):
	'''
	The probability density of t distribution.
	
	Parameters
	----------
	x : array_like
		Array of quantiles
	mu : mean
		The mean of normal distribution
	sig : sigma
		The standard deviation of normal distribution
	
	Returns
	-------
	logpdf : array_like
		Log of the normal probability density function evaluated at x
	
	'''
	return np.log(stats.t.pdf(x, loc = mu, scale = sig, df = len(x)-1))


#def dnorm(x, mu, sig):
#    return np.log(1/(sig * np.sqrt(2 * np.pi)) * np.exp(-(x - mu)**2 / (2 * sig**2)))

def dnorm(x, mu, sig):
	'''
	The probability density of normal distribution.
	
	Parameters
	----------
	x : array_like
		Array of quantiles
	mu : mean
		The mean of normal distribution
	sig : sigma
		The standard deviation of normal distribution
	
	Returns
	-------
	logpdf : array_like
		Log of the normal probability density function evaluated at x
	
	Notes
	------
	Do NOT set log_p to False in this file
		
	'''
	return np.log(stats.norm.pdf(x, mu, sig))

def dexp(x, l):
	'''
	The probability density of exponential distribution.
	
	Parameters
	----------
	x : array_like
		Array of quantiles
	l : lambda
		A common parameter for `expon`. such that ``pdf = lambda * exp(-lambda * x)``.
		The `scale` parameter in `expon` function corresponds to `scale = `1/lambda`
	log_p : bool, optional
		If True, return log(p-value)

	Returns
	-------
	logpdf : array_like
		Log of the exponential probability density function evaluated at x

	Notes
	------
	Do NOT set log_p to False in this file		
	'''
	
	if x > 0:
		return np.log(stats.expon.pdf(x, scale = 1/l))
	else:
		return 0

def like(s1, s2, para):
	'''
	Estimate the likelihood of observing data (s1 and s2) given parameters `para`

	Parameters
	----------
	s1 : array_like
		Beta values in group-1
	s2 : array_like
		Beta values in group-2
	para : list
		parameters consist of [mu1, sig1, mu2, sig2]
	
	Returns
	--------
	likelihood : float
		The log likelihood of observing data (s1 and s2) given parameters `para`
		
	Notes
	-----
	`log_p` in `dnorm` and `dexp` must be set to `True`
	'''
	[mu1, sig1, mu2, sig2] = para
	return np.sum(dt(s1, mu1, sig1)) + np.sum(dt(s2, mu2, sig2))

def prior(s1, s2, para):
	'''
	Probability of mean and std.
	
	Parameters
	----------
	s1 : array_like
		Beta values in group-1
	s2 : array_like
		Beta values in group-2
	para : list
		parameters consist of [mu1, sig1, mu2, sig2]
	
	Returns
	--------
	likelihood : float
		The log likelihood of `mean` which follows the normal distribution whose
		mean = pooled.mean, and whose std = pool.std * 1000. `Std` follows exponential
		distribution 
		
	Notes
	-----
	`log_p` in `dnorm` and `dexp` must be set to `True`		
	'''
	[mu1, sig1, mu2, sig2] = para
	pooled = np.append(s1, s2)
	prior_mean = pooled.mean()
	prior_std = 1000.0*pooled.std()
	
	return np.sum([dnorm(mu1, prior_mean, prior_std), dnorm(mu2, prior_mean, prior_std), dexp(sig1, 0.1), dexp(sig2, 0.1)])

def posterior(s1, s2, para):
	'''
	Parameters
	----------
	s1 : array_like
		Beta values in group-1
	s2 : array_like
		Beta values in group-2
	para : list
		parameters consist of [mu1, sig1, mu2, sig2]
	
	Returns
	-------
	likelihood : float
		The log likelihood of posterior distribution
			
	'''
	[mu1, sig1, mu2, sig2] = para
	return like(s1, s2, [mu1, sig1, mu2, sig2]) + prior(s1, s2, [mu1, sig1, mu2, sig2])

def computeHDI(chain, interval = .95):
	'''
	Compute 95% highest density interval (HDI)
	'''
	# sort chain using the first axis which is the chain
	chain.sort()
	# how many samples did you generate?
	nSample = chain.size    
	# how many samples must go in the HDI?
	nSampleCred = int(np.ceil(nSample * interval))
	# number of intervals to be compared
	nCI = nSample - nSampleCred
	# width of every proposed interval
	width = np.array([chain[i+nSampleCred] - chain[i] for  i in range(nCI)])
	# index of lower bound of shortest interval (which is the HDI) 
	best  = width.argmin()
	# put it in a dictionary
	#HDI   = {'Lower': chain[best], 'Upper': chain[best + nSampleCred], 'Width': width.min()}
	HDI = [chain[best], chain[best + nSampleCred]]
	return HDI
    

def beta_bayes(results, id, s1, s2, seed, niter = 10000, nburn_in = 500):
	'''
	https://stats.stackexchange.com/questions/130389/bayesian-equivalent-of-two-sample-t-test
	'''
	np.random.seed(seed)
	
	mu1_samples = []	#means sampled by MCMC for s1
	mu2_samples = []	#means sampled by MCMC for s2
	
	# run MCMC (Metropolis-Hastings's sampling algorithm)
	# Initialization: mu1, sig1, mu2, sig2
	parameters = np.array([np.mean(s1), np.std(s1), np.mean(s2), np.std(s2)])	
	increment = (s1.std() + s2.std())/10	#5 times smaller than the average std
	for iteration in np.arange(1,niter):
		candidate = parameters + np.random.normal(0, increment, 4)
		if candidate[1] < 0 or candidate[3] < 0:
			continue
		ratio = np.exp(posterior(s1, s2, candidate) - posterior(s1, s2, parameters))
		if np.random.uniform() < ratio:
			parameters = candidate
		if iteration < nburn_in:
			continue
		mu1_samples.append(parameters[0])
		mu2_samples.append(parameters[2])
	
	# calculate estimated means			
	mu1_samples = np.array(mu1_samples)
	mu2_samples = np.array(mu2_samples)
	est_mu1 = mu1_samples.mean()	#estimated mu1
	est_mu2 = mu2_samples.mean()	#estimated mu2
	
	# calculate probability
	diff = (mu1_samples - mu2_samples)
	diff_median = np.median(diff)
	if diff_median < 0:
		prob = np.mean(diff < 0)
	elif diff_median > 0:
		prob = np.mean(diff > 0)
	
	# calculate HDI
	diff_HDI_h, diff_HDI_l = computeHDI(diff)
	
	# CpG_ID, mean of group1, mean of group2, diff of mean, 95%HDI_low, HDI_high, probability
	results.append( [id, est_mu1, est_mu2, est_mu1 - est_mu2, diff_HDI_h, diff_HDI_l,  prob])

def beta_bayes_new(results, id, s1, s2, seed, niter = 10000, nburn_in = 500):
	'''
	https://stats.stackexchange.com/questions/130389/bayesian-equivalent-of-two-sample-t-test
	'''
	np.random.seed(seed)
	
	mu1_samples = []	#means sampled by MCMC for s1
	mu2_samples = []	#means sampled by MCMC for s2
	
	# run MCMC (Metropolis-Hastings's sampling algorithm)
	# Initialization: mu1, sig1, mu2, sig2
	parameters = np.array([np.mean(s1), np.std(s1), np.mean(s2), np.std(s2)])	
	increment = (s1.std() + s2.std())/10	#5 times smaller than the average std
	for iteration in np.arange(1,niter):
		candidate = parameters + np.random.normal(0, increment, 4)
		if candidate[1] < 0 or candidate[3] < 0:
			continue
		ratio = np.exp(posterior(s1, s2, candidate) - posterior(s1, s2, parameters))
		if np.random.uniform() < ratio:
			parameters = candidate
		if iteration < nburn_in:
			continue
		mu1_samples.append(parameters[0])
		mu2_samples.append(parameters[2])
	
	# calculate estimated means			
	mu1_samples = np.array(mu1_samples)
	mu2_samples = np.array(mu2_samples)
	est_mu1 = mu1_samples.mean()	#estimated mu1
	est_mu2 = mu2_samples.mean()	#estimated mu2
	
	# calculate probability
	diff = (mu1_samples - mu2_samples)
	diff_median = np.median(diff)
	if diff_median < 0:
		prob = np.mean(diff < 0)
	elif diff_median > 0:
		prob = np.mean(diff > 0)
	
	# calculate HDI
	diff_HDI_h, diff_HDI_l = computeHDI(diff)
	
	# CpG_ID, mean of group1, mean of group2, diff of mean, 95%HDI_low, HDI_high, probability
	results.append( [id, est_mu1, est_mu2, est_mu1 - est_mu2, diff_HDI_h, diff_HDI_l,  prob])

	
def test():
	np.random.seed(99)
	sample1 = np.random.normal(100, 10, 8)
	sample2 = np.random.normal(150, 15, 10)
	print (','.join([str(i) for i in sample1]))
	print (','.join([str(i) for i in sample2]))
	out = beta_bayes('test', sample1, sample2)
	print (out)
	
if __name__=='__main__':
	#test()
	
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing beta values with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). Except for the 1st row and 1st column, any non-numerical values will be considered as \"missing values\" and ignored. This file can be a regular text file or compressed file (*.gz, *.bz2) or accessible url.")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological group of each sample. It is a comma-separated 2 columns file with the 1st column containing sample IDs, and the 2nd column containing group IDs.  It must have a header row. Sample IDs should match to the \"Data file\". Note: Only for two group comparison.")
	parser.add_option("-n","--niter",action="store",type="int", default=5000,dest="n_iter",help="Iteration times when using MCMC Metropolis-Hastings's agorithm to draw samples from the posterior distribution. default=%default")
	parser.add_option("-b","--burnin",action="store",type="int", default=500,dest="n_burn",help="Number of samples to discard. Thes initial samples are usually not completely valid because the Markov Chain has not stabilized to the stationary distributio. default=%default.")
	parser.add_option("-p","--processor",action="store",type="int",dest="n_process",default=1,help="Number of processes. default=%default")
	parser.add_option("-s","--seed",action="store",type='int', dest="seed",default=99, help="The seed used by the random number generator. default=%default")
	parser.add_option("-o","--output",action="store",type="string", dest="out_file",help="Prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()
	#print (options.paired)
	#print (options.welch_ttest)
	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)

	if not (options.group_file):
		print (__doc__)
		parser.print_help()
		sys.exit(102)
				
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)	
	if options.n_iter <= 0:
		print ("--niter must be a positive integer")
		parser.print_help()		
		sys.exit(0)	
	if options.n_burn <= 0:
		print ("--burnin must be a positive integer")
		parser.print_help()	
		sys.exit(0)		
	if options.n_process <= 0:
		print ("--processor must be a positive integer")
		parser.print_help()		
		sys.exit(0)	
	
	np.random.seed(options.seed)
	printlog("Read group file \"%s\" ..." % (options.group_file))
	(ss,gs) = read_grp_file1(options.group_file)
	
	s2g = {}
	for s,g in zip(ss,gs):
		s2g[s] = g
	
	g2s = collections.defaultdict(list)
	for s,g in zip(ss, gs):
		g2s[g].append(s)
	
	group_IDs = sorted(g2s.keys())
	for g in group_IDs:
		print ("\tGroup %s has %d samples:" % (g, len(g2s[g])), file=sys.stderr)
		print ('\t\t' + ','.join(g2s[g]), file=sys.stderr)
	
	if len(group_IDs) != 2:
		printlog("You must have two groups!", file=sys.stderr)
		sys.exit(1)
	

	manager = Manager()
	results = manager.list()	 #list of list. shared variable between main() and beta_bayes(). #ID, group1.mean, group2.mean, prob
        	
	printlog("Read data file \"%s\" ..." % (options.input_file))
	line_num = 0
	p_count = 0
	jobs = []
	for l in ireader.reader(options.input_file):
		line_num += 1
		f = l.split()
		if line_num == 1:
			sample_IDs = f[1:]
			# check if sample ID matches
			for s in s2g:
				if s not in sample_IDs:
					printlog("Cannot find sample ID \"%s\" from file \"%s\"" % (s, options.input_file))
					sys.exit(3)
			g_IDs = [s2g[i] for i in sample_IDs]
			
		else:
			probe_ID = f[0]
			p_count += 1
			group1 = []	#beta values in group1
			group2 = [] #beta values in group2
			
			beta_values = f[1:]
			for g,b in zip(g_IDs, beta_values):
				#deal with non-numerical values
				try:
					b = float(b)
				except:
					continue
				if g == group_IDs[0]:
					group1.append(b)
				elif  g == group_IDs[1]:
					group2.append(b)
			
			group1 = np.array(group1)
			group2 = np.array(group2)			
			job_name = probe_ID
			p = Process(name = job_name,target = beta_bayes, args = (results, probe_ID, group1, group2, options.seed, options.n_iter, options.n_burn))
			p.start()
			jobs.append(p)
			
			if p_count == options.n_process:
				for proc in jobs: proc.join()	#tell the process to complete
				p_count = 0
				jobs = []
			print("Finish %d\r" % (line_num - 1),end = '', file=sys.stderr)
	for proc in jobs: proc.join()
	
	OUT = open(options.out_file + '.bayes.tsv','w')
	print ("\t".join(["ID", "mu1", "mu2", "mu_diff", "mu_diff (95% HDI)", "Probability"]), file=OUT)
	for r in results:
		print ("%s\t%f\t%f\t%f\t(%f,%f)\t%f" % (r[0], r[1],r[2],r[3], r[4], r[5],r[6]), file = OUT)
	OUT.close()
	
