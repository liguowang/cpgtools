dmc_Bayes.py
=============

Description
-----------

Different from statistical testing, this program tries to estimates "how different the
means between the two groups are" using Bayesian approach. An `MCMC <https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo>`_
is used to estimate the "means", "difference of means", "95% HDI (highest posterior density interval)",
and the posterior probability that the HDI does NOT include "0".

It is similar to John Kruschke's `BEST algorithm <http://www.indiana.edu/~kruschke/BEST/>`_
(Bayesian Estimation Supersedes T test)

**Notes**

- This program is much slower than T test due to MCMC (Markov chain Monte Carlo) step. 
  Running it with multiple threads is highly recommended.


Options
----------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file-INPUT_FILE
                        Data file containing beta values with the 1st row
                        containing sample IDs (must be unique) and the 1st
                        column containing CpG positions or probe IDs (must be
                        unique). Except for the 1st row and 1st column, any
                        non-numerical values will be considered as "missing
                        values" and ignored. This file can be a regular text
                        file or compressed file (.gz, .bz2).
  -g GROUP_FILE, --group-GROUP_FILE
                        Group file defining the biological group of each
                        sample. It is a comma-separated 2 columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs.  It must have a header
                        row. Sample IDs should match to the "Data file". Note:
                        Only for two group comparison.
  -n N_ITER, --niter-N_ITER
                        Iteration times when using MCMC Metropolis-Hastings's
                        agorithm to draw samples from the posterior
                        distribution. default-5000
  -b N_BURN, --burnin-N_BURN
                        Number of samples to discard. Thes initial samples are
                        usually not completely valid because the Markov Chain
                        has not stabilized to the stationary distributio.
                        default-500.
  -p N_PROCESS, --processor-N_PROCESS
                        Number of processes. default-1
  -s SEED, --seed-SEED  The seed used by the random number generator.
                        default-99
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.



Input files (examples)
------------------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
- `test_05_TwoGroup.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_

Command
---------

::

 $  dmc_Bayes.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -p 10 -o dmc_output                        

Output files
-----------------

- **dmc_output.bayes.tsv**: this file consists of 6 columns:
 
 1. ID : CpG ID
 2. *mu1* : Mean methylation level estimated from group1
 3. *mu2* : Mean methylation level estimated from gropu2
 4. *mu_diff* : Difference between mu1 and mu2
 5. *mu_diff* (95% HDI) : 95% of "High Density Interval" of *mu_diff*. The HDI indicates which
    points of a distribution are most credible. This interval spans 95% of *mu_diff*'s
    distribution. 
 6. The probability that *mu1* and *mu2* are different. 
    
::

 $head -10 dmc_output.bayes.tsv
 
 ID	mu1	mu2	mu_diff	mu_diff (95% HDI)	Probability
 cg00001099	0.775209	0.795404	-0.020196	(-0.065148,0.023974)	0.811024
 cg00000363	0.610565	0.469523	0.141042	(0.030769,0.232965)	0.994665
 cg00000884	0.845973	0.873761	-0.027787	(-0.051976,-0.004398)	0.984882
 cg00000714	0.190868	0.199233	-0.008365	(-0.030071,0.014006)	0.816141
 cg00000957	0.772905	0.827528	-0.054623	(-0.092116,-0.016465)	0.995327
 cg00000292	0.748394	0.766326	-0.017932	(-0.051286,0.012583)	0.889729
 cg00000807	0.729162	0.683732	0.045430	(-0.001523,0.086588)	0.981551
 cg00000721	0.935903	0.935080	0.000823	(-0.013210,0.018628)	0.508686
 cg00000948	0.898609	0.897536	0.001073	(-0.020663,0.026813)	0.518238

