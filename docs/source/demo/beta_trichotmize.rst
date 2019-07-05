beta_trichotmize.py
====================

Description
--------------
Rather than using hard threshold to call "methylated" or "unmethylated" CpGs or regions, 
this program uses probability approach (Bayesian Gaussian Mixture model) to trichotmize
beta values into three status:

Un-methylated : labeled as "0" in result file
	Both the homologous chromosomes (i.e. The maternal and paternal chromosomes) are unmethylated. 
Semi-methylated : labeled as "1" in result file
	Only one of the homologous chromosomes is methylated. This is also called allele-specific
	methylation or "imprinting". Note: this is different from **hemimethylation**, which refers
	to "one of two (complementary) strands is methylated".  
Full-methylated : labeled as "2" in result file
	Both the homologous chromosomes (i.e. The maternal and paternal chromosomes) are methylated. 
unassigned : labeled as "-1" in result file
	CpGs failed to assigned to the three categories above.
	
Algorithm
---------
As described above, in somatic cells, most CpGs can be grouped into 3 categories including
"Un-methylated", "Semi-methylated (imprinted)" and "Full-methylated". Therefore, the
Beta distribution of CpGs can be considered as the mixture of 3 Gaussian distributions
(i.e. components). **beta_trichotmize.py** first estimates the parameters (mu1, mu2, mu3)
and (s1, s2, s3) of the 3 components using expectation–maximization (EM) algorithm, then it 
calculates the posterior probabilities ( *p0*, *p1*, and *p2*) of each component given
the beta value of a CpG. 


*p0*
	the probability that the CpG belongs to **un-methylated** component. 
*p1*
	the probability that the CpG belongs to **semi-methylated**  component. 
*p2*
	the probability that the CpG belongs to **full-methylated** component. 

The classification will be made using rules:

::

 if p0 == max(p0, p1, p2):
 	un-methylated
 elif p2 == max(p0, p1, p2):
 	full-methylated
 elif p1 == max(p0, p1, p2):
 	if p1 >= prob_cutoff:
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
