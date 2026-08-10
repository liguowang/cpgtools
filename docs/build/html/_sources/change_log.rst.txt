CpGtools Release History
========================

Version 3.0.0
-------------

Major release introducing a comprehensive DNA methylation missing-value
imputation framework.

New features

* Added :code:`beta_impute.py`, a unified command-line framework for DNA
  methylation missing-value analysis.
* Added multiple imputation algorithms, including:

  * Constant replacement
  * Mean, median, minimum, and maximum imputation
  * Random-value imputation
  * Moving-window imputation
  * K-nearest neighbors (KNN)
  * Reference-based KNN
  * Iterative Buck regression
  * Iterative Random Forest
  * SoftImpute matrix completion
  * MOREL block-wise imputation
  * Genomic nearest-neighbor (GNN) imputation

* Added utilities for

  * generating synthetic methylation matrices,
  * inserting missing values into existing matrices,
  * summarizing missing values, and
  * evaluating imputation accuracy against a truth matrix using MAE, RMSE,
    and R².

Documentation

* Added comprehensive documentation for :code:`beta_impute.py`,
  including installation instructions, algorithm selection guidance,
  and usage examples.

Version 2.0.4
-------------

New features

* Added :code:`beta_combat.py`.

Version 2.0.3
-------------

Bug fixes

* Fixed an issue in the ANOVA workflow where p-values and adjusted p-values
  were reported as missing for all CpGs.

Version 2.0.1
-------------

New features

* Added :code:`predict_sex.py`.
* Added :code:`beta_selectNBest.py`.

Version 1.10.0
--------------

New features

* Added :code:`beta_UMAP.py`.

Version 1.0.8
-------------

Bug fixes

* Fixed an issue in :code:`beta_tSNE.py` and :code:`beta_PCA.py` when sample
  identifiers were numeric.

Version 1.0.7
-------------

New features

* Added :code:`CpG_density_gene_centered.py`.

Version 1.0.2
-------------

New features

* Added :code:`beta_tSNE.py` for t-distributed stochastic neighbor embedding
  (t-SNE) analysis of DNA methylation samples.

Version 1.0.1
-------------

New features

* Added :code:`CpG_anno_position.py` for annotating CpGs using pre-built or
  user-supplied genomic annotation files.

Version 1.0.0
-------------

Initial public release.
