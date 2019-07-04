beta_trichotmize.py
====================

Description
--------------
Rather than using hard threshold to call "methylated" or "unmethylated" CpGs or regions, 
this program uses probability approach (Bayesian Gaussian Mixture model) to trichotmize
beta values into three status:

- Un-methylated (labeled as "0" in result file)
- Semi-methylated (labeled as "1" in result file)
- Full-methylated (labeled as "2" in result file)
- unassigned (labeled as "-1" in result file)

Basically, GMM will first calculate probability *p0*, *p1*, and *p2* for each CpG based
on its beta value:

*p0*
	the probability that the CpG is un-methylated
*p1*
	the probability that the CpG is semi-methylated
*p2*
	the probability that the CpG is full-methylated

The classification will be made using rules:

::

 if p0 -- max(p0, p1, p2):
 	un-methylated
 elif p2 -- max(p0, p1, p2):
 	full-methylated
 elif p1 -- max(p0, p1, p2):
 	if p1 >- prob_cutoff:
 		semi-methylated
 	else:
 	 	unknown/unassigned

Input files (examples)
------------------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_

Command
--------
::

 $beta_trichotmize.py -i test_05_TwoGroup.tsv -r

Below histogram and piechart showed the proportion of CpGs assigned to "Un-methylated", "Semi-methylated" and "Full-methylated". 

.. image:: ../_static/trichotmize.png
   :height: 650 px
   :width: 650 px
   :scale: 100 %  
