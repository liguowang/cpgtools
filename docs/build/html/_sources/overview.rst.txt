Overview
========

**CpGtools** is a Python package providing a comprehensive collection of
command-line tools for the analysis of DNA methylation data. The package
supports methylation data generated from both microarray and sequencing
platforms, including:

* `Illumina HumanMethylation450 BeadChip (450K) <https://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit.html>`_
* `Illumina MethylationEPIC BeadChip (EPIC, 850K) <https://www.illumina.com/products/by-type/microarray-kits/infinium-methylationepic-beadchip-kit.html>`_
* `Illumina Infinium MethylationEPIC v2.0 (EPIC v2, 930K) <https://www.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html>`_
* `Reduced Representation Bisulfite Sequencing (RRBS) <https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/rrbs-seq-scrrbs.html>`_
* Whole-Genome Bisulfite Sequencing (WGBS)

CpGtools provides utilities for

* annotating CpGs and genomic regions,
* summarizing and visualizing DNA methylation profiles,
* identifying differentially methylated CpGs,
* imputing missing values using a wide range of statistical and
  machine-learning algorithms, and
* predicting biological phenotypes from DNA methylation profiles.

The modules are organized into the following functional categories:

#. **CpG Position Analysis** – annotate CpGs and characterize their genomic
   distribution.
#. **CpG Signal Analysis** – summarize, normalize, transform, and visualize
   DNA methylation signals.
#. **Missing-value Imputation** – simulate, evaluate, and impute missing
   values in DNA methylation matrices.
#. **Differential Methylation Analysis** – identify differentially methylated
   CpGs using statistical and regression-based approaches.
#. **Predictive Analysis** – predict biological phenotypes and sample
   characteristics from DNA methylation profiles.

Most tools are designed for command-line use and follow a consistent interface.
Detailed documentation, usage examples, and command-line options are available
on the corresponding documentation pages.

CpG Position Analysis
---------------------

These modules annotate CpGs, summarize their genomic distribution, and
characterize their relationships with genes and other genomic features.

.. list-table:: CpG position analysis modules
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Description

   * - ``CpG_aggregation.py``
     - Aggregate CpG methylation values within user-defined genomic regions
       (e.g., promoters, CpG islands, gene bodies, and exons).

   * - ``CpG_anno_position.py``
     - Annotate CpGs according to their genomic coordinates using built-in
       or user-supplied annotation files.

   * - ``CpG_anno_probe.py``
     - Annotate Illumina HumanMethylation450 (450K), EPIC (850K), and
       EPIC v2 probes with genomic and functional information.

   * - ``CpG_density_gene_centered.py``
     - Calculate and visualize CpG density across gene-centered genomic
       regions, including upstream, gene body, and downstream regions.

   * - ``CpG_distrb_chrom.py``
     - Summarize the chromosomal distribution of CpGs.

   * - ``CpG_distrb_gene_centered.py``
     - Summarize the distribution of CpGs across gene-centered genomic
       features.

   * - ``CpG_distrb_region.py``
     - Summarize the distribution of CpGs across user-defined genomic
       regions.

   * - ``CpG_logo.py``
     - Generate DNA sequence logos and position weight matrices for CpGs
       and their flanking sequences.

   * - ``CpG_to_gene.py``
     - Assign CpGs to putative target genes using an algorithm similar to
       `GREAT <http://great.stanford.edu/public/html/>`_.

CpG Signal Analysis
-------------------

These modules summarize, transform, normalize, visualize, and explore DNA
methylation signals across CpGs and samples. Most modules support both
DNA methylation **beta-value** and **M-value** matrices.

.. list-table:: CpG signal analysis modules
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Description

   * - ``beta_combat.py``
     - Correct batch effects using the
       `ComBat <https://pubmed.ncbi.nlm.nih.gov/16632515/>`_ algorithm.

   * - ``beta_jitter_plot.py``
     - Generate jitter plots (strip charts) and bean plots to visualize the
       distribution of methylation values for each sample.

   * - ``beta_m_conversion.py``
     - Convert DNA methylation beta values to M values, or vice versa.

   * - ``beta_PCA.py``
     - Perform principal component analysis (PCA) to explore sample
       relationships and identify major sources of variation.

   * - ``beta_profile_gene_centered.py``
     - Calculate and visualize average methylation profiles across
       gene-centered genomic regions.

   * - ``beta_profile_region.py``
     - Calculate and visualize average methylation profiles across
       user-defined genomic regions.

   * - ``beta_selectNBest.py``
     - Select the most informative CpGs using ANOVA, mutual information,
       or chi-square feature selection.

   * - ``beta_stacked_barplot.py``
     - Generate stacked bar plots showing the proportions of CpGs within
       different methylation intervals for each sample.

   * - ``beta_stats.py``
     - Summarize methylation statistics for CpGs within genomic regions.

   * - ``beta_topN.py``
     - Select the most variable CpGs according to standard deviation.

   * - ``beta_trichotomize.py``
     - Classify beta values into unmethylated, partially methylated,
       methylated, or unassigned states using a Bayesian Gaussian mixture
       model.

   * - ``beta_tSNE.py``
     - Perform t-distributed stochastic neighbor embedding (t-SNE) for
       visualization of sample relationships.

   * - ``beta_UMAP.py``
     - Perform Uniform Manifold Approximation and Projection (UMAP) for
       visualization of sample relationships.

Differential Methylation Analysis
---------------------------------

These modules identify **differentially methylated CpGs (DMCs)** between
experimental or biological groups using statistical tests, regression models,
and Bayesian methods. Different modules are available for Illumina methylation
arrays (450K, EPIC, and EPIC v2) and bisulfite sequencing data (RRBS/WGBS).

.. list-table:: Differential methylation analysis modules
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Description

   * - ``dmc_Bayes.py``
     - Identify differentially methylated CpGs using a Bayesian statistical
       model. Designed for Illumina methylation array data.

   * - ``dmc_bb.py``
     - Identify differentially methylated CpGs from RRBS/WGBS count data
       using a beta-binomial model.

   * - ``dmc_fisher.py``
     - Identify differentially methylated CpGs from RRBS/WGBS count data
       using Fisher's exact test.

   * - ``dmc_glm.py``
     - Identify differentially methylated CpGs using a generalized linear
       model (GLM). Designed for Illumina methylation array data.

   * - ``dmc_logit.py``
     - Identify differentially methylated CpGs from RRBS/WGBS count data
       using logistic regression.

   * - ``dmc_nonparametric.py``
     - Identify differentially methylated CpGs using nonparametric tests,
       including the Mann–Whitney U test for two-group comparisons and the
       Kruskal–Wallis test for multiple-group comparisons.

   * - ``dmc_ttest.py``
     - Identify differentially methylated CpGs using Student's *t*-test.
       Designed for Illumina methylation array data.

Predictive Analysis
-------------------

These modules predict biological phenotypes or sample characteristics from
DNA methylation profiles. Additional predictive models will be incorporated
in future releases.

.. list-table:: Predictive analysis modules
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Description

   * - ``predict_sex.py``
     - Predict the biological sex of a sample using DNA methylation patterns
       on sex chromosomes and CpGs associated with genomic imprinting.

Missing-value Imputation
------------------------

These modules simulate, summarize, evaluate, and impute missing values in DNA
methylation matrices. The framework supports both **beta-value** and
**M-value** matrices and provides a broad range of statistical,
machine-learning, and genomics-informed imputation algorithms.

.. list-table:: Missing-value imputation modules
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Description

   * - ``beta_impute.py``
     - A comprehensive framework for DNA methylation missing-value analysis.
       Supports data simulation, missing-value insertion, multiple
       imputation algorithms, and quantitative evaluation of imputation
       accuracy.

The :code:`beta_impute.py` command provides the following functionality:

* **Simulation and benchmarking**

  * Generate synthetic DNA methylation matrices with user-defined dimensions
    and missing-value rates.
  * Insert missing values into existing matrices for benchmarking and method
    comparison.

* **Missing-value assessment**

  * Count missing values for each CpG and sample.
  * Remove incomplete rows or columns when required.

* **Missing-value imputation**

  * Simple statistical methods (constant, mean, median, minimum, maximum,
    and random-value imputation).
  * Moving-window imputation.
  * K-nearest neighbors (KNN).
  * Reference-based KNN using an external complete methylation matrix.
  * Iterative Buck regression.
  * Iterative Random Forest regression.
  * SoftImpute low-rank matrix completion.
  * MOREL for systematic block-wise missingness.
  * Genomic nearest-neighbor (GNN) imputation using neighboring CpGs.

* **Performance evaluation**

  * Compare imputed values with a complete reference (truth) matrix.
  * Report Mean Absolute Error (MAE), Root Mean Squared Error (RMSE), and
    coefficient of determination (R²) for the imputed values.

For detailed documentation, available algorithms, command-line options, and
examples, see :doc:`demo/beta_impute`.
