Overview
=========

CpGtools package provides a number of Python programs to annotate, QC, visualize, and
analyze DNA methylation data generated from Illumina
`HumanMethylation450 BeadChip (450K) <https://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit.html>`_ /
`MethylationEPIC BeadChip (850K) <https://www.illumina.com/documents/products/datasheets/datasheet_CytoSNP850K_POP.pdf>`_ array or
`RRBS / WGBS <https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/rrbs-seq-scrrbs.html>`_.

These programs can be divided into three groups:

- CpG position analysis modules
- CpG signal analysis modules
- Differential CpG analysis modules

CpG position analysis modules
-----------------------------
These modules are primarily used to analyze CpG's genomic locations. 

+------------------------------+-------------------------------------------------------------------+
|Name                          |Description                                                        |
+------------------------------+-------------------------------------------------------------------+
|CpG_anno_probe.py             |add comprehensive annotation information to each 450K/850K probe   |
|                              |ID.                                                                |
+------------------------------+-------------------------------------------------------------------+
|CpG_aggregation.py            |Aggregate proportion values of a list of CpGs that located in give |
|                              |genomic regions.                                                   |
+------------------------------+-------------------------------------------------------------------+
|CpG_distrb_chrom.py           |Calculates the distribution of CpG over chromosomes.               |
+------------------------------+-------------------------------------------------------------------+
|CpG_distrb_gene_centered.py   |Calculates the distribution of CpG over gene-centered genomic      |
|                              |regions including 'Coding exons', 'UTR exons', 'Introns',          |
|                              |'Upstream intergenic regions', and 'Downsteam intergenic regions'. |
+------------------------------+-------------------------------------------------------------------+
|CpG_distrb_region.py          |Calculates the distribution of CpG over user-specified genomic     |
|                              |regions (such as promoters, enhancers).                            |
+------------------------------+-------------------------------------------------------------------+
|CpG_logo.py                   |Generates DNA motif logo for a given set of CpGs (to visualize     |
|                              |the genomic context of    these CpGs).                             |
+------------------------------+-------------------------------------------------------------------+
|CpG_to_gene.py                |Assigns CpGs to their putative target genes. Follows the "Basel    |
|                              |plus extension" rules used by the `GREAT <http://great.stanford.edu|
|                              |/public/html/index.php>`_ algorithm.                               |
+------------------------------+-------------------------------------------------------------------+

CpG signal analysis modules
----------------------------
These modules are primarily used to analyze CpG's DNA methylation beta values 

+------------------------------+-------------------------------------------------------------------+
|Name                          |Description                                                        |
+------------------------------+-------------------------------------------------------------------+
|beta_PCA.py                   |Performs `PCA <https://en.wikipedia.org/wiki/Principal_component_  |
|                              |analysis>`_ for samples.                                           |
+------------------------------+-------------------------------------------------------------------+
|beta_jitter_plot.py           |Generates jitter plot and bean plot for each sample (column).      |
+------------------------------+-------------------------------------------------------------------+
|beta_m_conversion.py          |Converts Beta-value into M-value or vice versa.                    |
+------------------------------+-------------------------------------------------------------------+
|beta_profile_gene_centered.py |Calculates the methylation profile (i.e. average beta value) for   |
|                              |gene-centered genomic regions.                                     |
+------------------------------+-------------------------------------------------------------------+
|beta_profile_region.py        |Calculates the methylation profile (i.e. average beta value) for   |
|                              |user specified genomic regions.                                    |
+------------------------------+-------------------------------------------------------------------+
|beta_stacked_barplot.py       |Creates stacked barplot for each sample. The stacked barplot       |
|                              |showing the proportions of CpGs whose beta values are falling into |
|                              |these 4 ranges: [0.00,  0.25], [0.25,  0.50], [0.50,  0.75], and   |
|                              |[0.75,  1.00].                                                     |
+------------------------------+-------------------------------------------------------------------+
|beta_m_conversion.py          |Gives basic information of CpGs located in genomic regions. These  |
|                              |information include "Number of CpGs", "Min methylation level",     |
|                              |"Max methylation level", "Mean methylation level across all CpGs", |
|                              |"Median methylation level across all CpGs" and "Standard deviation"|
+------------------------------+-------------------------------------------------------------------+
|beta_m_conversion.py          |This program picks the N most variable CpGs from the input file.   |
|                              |The result file can be used for PCA or clustering analysis.        |
+------------------------------+-------------------------------------------------------------------+
|beta_trichotmize.py           |This program uses `Bayesian Gaussian Mixture model <https://scikit-|
|                              |learn.org/stable/modules/generated/sklearn.mixture.BayesianGaussian|
|                              |Mixture.html>`_ to trichotmize beta values into three status:      |
|                              |"Un-methylated","Semi-methylated", "Full-methylated", and          |
|                              |"unassigned"                                                       |
+------------------------------+-------------------------------------------------------------------+

Differential CpG analysis modules
----------------------------------
These modules are primarily used to identify differentially methylated CpGs

+------------------------------+-------------------------------------------------------------------+
|Name                          |Description                                                        |
+------------------------------+-------------------------------------------------------------------+
|dmc_Bayes.py                  |Different from statistical testing, this program tries to estimates|
|                              |"how different the means between the two groups are" using Bayesian|
|                              |approach. An `MCMC <https://en.wikipedia.org/wiki/Markov_chain_    |
|                              |Monte_Carlo>`_ is used to estimate the "means", "difference of     |
|                              |means", "95% HDI (highest posterior density interval)", and the    |
|                              |posterior probability that the HDI does NOT include "0". It is     |
|                              |similar to John Kruschke's `BEST <(http://www.indiana.edu/~kruschke|
|                              |/BEST/)>`_ algorithm.                                              |
+------------------------------+-------------------------------------------------------------------+
|dmc_bb.py                     |This program performs differential CpG analysis using the `beta    |
|                              |binomial model <https://en.wikipedia.org/wiki/Beta-binomial        |
|                              |_distribution>`_ based on methylation proportions (in the form of  |
|                              |"c,n", where "c" indicates "number of reads with methylated C" and |
|                              |"n" indicates "Number of total reads".                             |
+------------------------------+-------------------------------------------------------------------+
|dmc_fisher.py                 |This program performs differential CpG analysis using Fisher's     |
|                              |Exact Test. It only applies to two sample comparison with no       |
|                              |replicates.                                                        |
+------------------------------+-------------------------------------------------------------------+
|dmc_glm.py                    |This program performs differential CpG analysis using the linear   |
|                              |regression model.                                                  |
+------------------------------+-------------------------------------------------------------------+
|dmc_logit.py                  |This program performs differential CpG analysis using the logistic |
|                              |regression model.                                                  |
+------------------------------+-------------------------------------------------------------------+
|dmc_nonparametric.py          |This program performs differential CpG analysis using `Mann-Whitney|
|                              |test <https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test>`_ |
|                              |for two group comparison, and the `Kruskal-Wallis H-test <https:// |
|                              |en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_  |
|                              |variance>`_ for multiple groups comparison.                        |
+------------------------------+-------------------------------------------------------------------+
|dmc_ttest.py                  |This program performs differential CpG analysis using the T test   |
|                              |for two group comparison, and `ANOVA <https://en.wikipedia.org/    |
|                              |wiki/Analysis_of_variance>`_ for multiple groups comparison.       |
+------------------------------+-------------------------------------------------------------------+
