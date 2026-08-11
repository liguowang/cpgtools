CpGtools Release History
========================

Version 3.0.0
-------------

Version 3.0.0 is a major release that expands CpGtools with epigenetic-clock
analysis, missing-value imputation, cell-type deconvolution, and a modernized
command-line and packaging structure.

New features
~~~~~~~~~~~~

* Added ``epical``, a unified command-line interface for calculating DNA
  methylation age and related aging measures using multiple published
  epigenetic clocks.
* Added ``beta_impute``, a unified framework for DNA methylation
  missing-value analysis and imputation.
* Added ``beta_deconvolution`` for estimating cell-type proportions from
  DNA methylation profiles.
* Added imputation methods including:

  * constant, mean, median, minimum, maximum, and random-value imputation;
  * moving-window and K-nearest-neighbor (KNN) imputation;
  * reference-based KNN;
  * iterative Buck regression and Random Forest regression;
  * SoftImpute matrix completion;
  * MOREL block-wise imputation; and
  * genomic nearest-neighbor (GNN) imputation.

* Added utilities for generating synthetic methylation matrices, inserting
  missing values, summarizing missingness, and evaluating imputation accuracy
  against a truth matrix using MAE, RMSE, and :math:`R^2`.

Command-line and packaging changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Standardized command-line programs as installed console commands.
* Reorganized command-line modules under ``cpgmodule.cli``.
* Refactored ``epical`` into a registry-based CLI architecture.
* Added automated tests for ``epical`` and ``beta_impute``.
* Modernized package configuration using ``pyproject.toml`` and the ``src``
  package layout.
* Updated dependency and license metadata for modern Python packaging.

Documentation
~~~~~~~~~~~~~

* Expanded documentation for ``epical``, ``beta_impute``, and
  ``beta_deconvolution``.
* Clarified input-file conventions, including the default methylation-matrix
  orientation of CpGs in rows and samples in columns.
* Updated installation and command-line usage documentation.


Version 2.0.4
-------------

New features
~~~~~~~~~~~~

* Added ``beta_combat`` for batch-effect correction.


Version 2.0.3
-------------

Bug fixes
~~~~~~~~~

* Fixed an issue in the ANOVA workflow where p-values and adjusted p-values
  were reported as missing for all CpGs.


Version 2.0.1
-------------

New features
~~~~~~~~~~~~

* Added ``predict_sex``.
* Added ``beta_selectNBest``.


Version 1.10.0
--------------

New features
~~~~~~~~~~~~

* Added ``beta_UMAP``.


Version 1.0.8
-------------

Bug fixes
~~~~~~~~~

* Fixed an issue in ``beta_tSNE`` and ``beta_PCA`` when sample identifiers
  were numeric.


Version 1.0.7
-------------

New features
~~~~~~~~~~~~

* Added ``CpG_density_gene_centered``.


Version 1.0.2
-------------

New features
~~~~~~~~~~~~

* Added ``beta_tSNE`` for t-distributed stochastic neighbor embedding
  (t-SNE) analysis of DNA methylation samples.


Version 1.0.1
-------------

New features
~~~~~~~~~~~~

* Added ``CpG_anno_position`` for annotating CpGs using pre-built or
  user-supplied genomic annotation files.


Version 1.0.0
-------------

Initial public release.
