beta_impute.py
==============

``beta_impute.py`` is a command-line framework for detecting, simulating, and
imputing missing values in DNA methylation matrices. Unlike most CpGtools
commands, which perform a single task, it provides multiple **subcommands**
for data simulation, missing-value assessment, imputation using a variety of
algorithms, and evaluation of imputation accuracy.

The program supports both DNA methylation beta-value and M-value matrices.
It is designed primarily for Illumina methylation arrays (e.g., 450K, EPIC,
and EPICv2), where **CpGs are represented as rows and samples as columns**,
but is equally applicable to reduced representation bisulfite sequencing
(RRBS), whole-genome bisulfite sequencing (WGBS), and other methylation
profiling platforms.

Highlights
----------

#. **Comprehensive collection of imputation methods.** Supports simple
   statistical methods (e.g., mean, median, moving window) as well as advanced
   machine-learning and matrix-completion approaches, including K-nearest
   neighbors (KNN), iterative regression (Buck), random forests, SoftImpute,
   deep neural networks (DNN), and genomic nearest-neighbor (GNN)
   imputation.

#. **Handles diverse missingness patterns.** Supports imputation of both
   sporadic (random) missing values and systematic (structural) missingness,
   such as large blocks of missing values caused by differences between
   methylation array platforms or sample groups.

#. **End-to-end benchmarking framework.** Provides utilities to generate
   simulated methylation matrices or insert missing values into existing
   matrices, impute missing values using multiple algorithms, and evaluate
   imputation accuracy against a complete reference matrix using metrics such
   as MAE, RMSE, and R². This enables users to compare methods and optimize
   algorithm parameters (e.g., tree depth, number of iterations, and the
   number of nearest neighbors) for their own data.


Command Summary
---------------

The available subcommands fall into three categories:

* **Simulation**
    * ``toy`` – generate a synthetic methylation matrix.
    * ``insertna`` – insert missing values into an existing matrix.

* **Utilities**
    * ``countna`` – summarize missing values by CpG and sample.
    * ``dropna`` – remove incomplete rows or columns.

* **Imputation**
    * ``constant``
    * ``mean``
    * ``median``
    * ``min``
    * ``max``
    * ``rand``
    * ``mw``
    * ``knn``
    * ``refknn``
    * ``buck``
    * ``rf``
    * ``softimpute``
    * ``morel``
    * ``gnn``

Run

.. code-block:: bash

   beta_impute.py <subcommand> --help

to display the command-line options for a specific subcommand.

Available Imputation Methods
-----------------------------

The table below summarizes the recommended use cases for each method.


.. list-table:: Implemented imputation methods
   :header-rows: 1
   :widths: 12 18 18 18 34

   * - Method
     - Recommended use
     - Matrix
     - Missingness
     - Notes

   * - ``constant``
     - Replace missing values with a fixed value.
     - Beta, M
     - Both
     - Useful for testing or special downstream workflows. Not recommended
       for biological analyses.

   * - ``mean``
     - Low levels of sporadic missing values.
     - Beta, M
     - Sporadic
     - Fast and simple baseline method.

   * - ``median``
     - Sporadic missing values with outliers.
     - Beta, M
     - Sporadic
     - More robust to extreme values than mean imputation.

   * - ``min``
     - Conservative replacement.
     - Beta, M
     - Sporadic
     - Mainly useful for sensitivity analyses.

   * - ``max``
     - Conservative replacement.
     - Beta, M
     - Sporadic
     - Rarely used in practice.

   * - ``rand``
     - Simulation or benchmarking.
     - Beta, M
     - Sporadic
     - Preserves the empirical distribution while introducing randomness.

   * - ``mw``
     - Missing values surrounded by nearby observations.
     - Beta, M
     - Sporadic
     - Assumes neighboring CpGs (or samples) have similar values.

   * - ``knn``
     - General sporadic missing values.
     - Beta, M
     - Sporadic
     - Good default choice for most DNA methylation data sets.

   * - ``refknn``
     - External reference matrix available.
     - Beta, M
     - Sporadic
     - Often improves accuracy when a large, biologically similar reference
       cohort is available.

   * - ``buck``
     - Moderate missingness with approximately linear relationships.
     - Beta, M
     - Sporadic
     - Iterative linear regression. May predict values outside the beta-value
       range.

   * - ``rf``
     - Complex nonlinear relationships.
     - Beta, M
     - Sporadic
     - Usually more accurate than linear regression but computationally more
       expensive. For beta-value matrices, consider specifying
       ``--min-value 0`` and ``--max-value 1``.

   * - ``softimpute``
     - Approximately low-rank matrices.
     - Beta, M
     - Sporadic
     - Matrix-completion method suitable for large data sets. For beta-value
       matrices, consider specifying ``--min-value 0`` and
       ``--max-value 1``.

   * - ``morel``
     - Systematic block-wise missingness between sample groups.
     - Beta, M
     - Block-wise
     - Designed for platform- or group-specific missingness. Supports both
       Random Forest and deep neural network models. For beta-value matrices,
       consider specifying ``--min-value 0`` and ``--max-value 1``.

   * - ``gnn``
     - CpGs with nearby genomic neighbors.
     - Beta, M
     - Sporadic
     - Uses neighboring CpGs based on genomic distance and, optionally,
       shared cis-regulatory elements (CREs).



**General recommendations**

* For **sporadic missing values**, start with ``knn``. It generally provides a
  good balance between accuracy and computational efficiency.

* When a **large, complete reference cohort** is available, consider
  ``refknn`` because biologically similar samples often improve imputation
  accuracy.

* For **systematic block-wise missingness** (for example, CpGs absent from one
  assay or one sample group), use ``morel`` rather than standard KNN-based
  methods.

* If **genomic annotations** are available and neighboring CpGs are expected
  to have correlated methylation levels, ``gnn`` can provide biologically
  meaningful estimates.

* ``rf`` generally models nonlinear relationships better than linear
  regression (``buck``), although it requires substantially more computation.

* ``softimpute`` is appropriate when the methylation matrix can be reasonably
  approximated by a low-rank structure.

* The statistical methods (``mean``, ``median``, ``min``, and ``max``) are
  computationally inexpensive and useful as baseline approaches for comparison,
  but they usually do not achieve the accuracy of model-based methods.


Utility Commands
-----------------

In addition to missing-value imputation, ``beta_impute.py`` provides several
utilities for data simulation, quality assessment, and preprocessing.

.. list-table:: Utility commands
   :header-rows: 1
   :widths: 18 32 50

   * - Command
     - Purpose
     - Description

   * - ``toy``
     - Generate a toy methylation matrix
     - Creates a synthetic DNA methylation matrix with user-specified
       dimensions and randomly inserted missing values. Useful for
       demonstrations, benchmarking, and software testing.

   * - ``insertna``
     - Simulate missing values
     - Randomly inserts missing values into an existing methylation matrix
       until a specified number of missing entries is reached. Useful for
       evaluating and comparing imputation methods.

   * - ``countna``
     - Summarize missing values
     - Reports the number of missing values for each CpG (row) and each
       sample (column), helping users assess the extent and distribution of
       missingness.

   * - ``dropna``
     - Remove incomplete data
     - Removes rows (CpGs) or columns (samples) containing one or more
       missing values to produce a complete matrix for downstream analyses.



Evaluating imputation accuracy
------------------------------

Most imputation commands support the optional ``--truth`` argument, which
specifies a complete reference matrix containing the true methylation values.
Only positions that were originally missing in the input matrix are evaluated.
Entries that remain missing after imputation are excluded automatically.

* **MAE (Mean Absolute Error)** – average absolute difference between the
  imputed and true values.
* **RMSE (Root Mean Squared Error)** – emphasizes larger imputation errors.
* **R² (Coefficient of Determination)** – measures the agreement between the
  imputed and true values.

Only positions that were originally missing in the input matrix are evaluated.
Entries that remain missing after imputation are excluded from the reported
metrics.

Examples
--------

The examples below illustrate common workflows. Additional options are
available for each subcommand; use

Generate a toy methylation matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a 1,000 × 20 DNA methylation matrix with 10% missing values.

.. code-block:: bash

   beta_impute.py toy toy.tsv \
       --rows 1000 \
       --cols 20 \
       --missingness 0.1

Count missing values
~~~~~~~~~~~~~~~~~~~~

Generate reports summarizing the number of missing values for each CpG and
each sample.

.. code-block:: bash

   beta_impute.py countna methylation.tsv

Insert missing values
~~~~~~~~~~~~~~~~~~~~~

Randomly insert missing values into an existing methylation matrix until the
matrix contains 5,000 missing entries.

.. code-block:: bash

   beta_impute.py insertna \
       methylation.tsv \
       methylation.missing.tsv \
       5000

Mean imputation
~~~~~~~~~~~~~~~

Replace missing values using the mean of each CpG (row).

.. code-block:: bash

   beta_impute.py mean \
       methylation.tsv \
       methylation.mean.tsv

K-nearest-neighbor imputation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Impute sporadic missing values by searching for similar samples.

.. code-block:: bash

   beta_impute.py knn \
       methylation.tsv \
       methylation.knn.tsv \
       --neighbors 5 \
       --weights distance

Reference-based KNN
~~~~~~~~~~~~~~~~~~~

Use a complete external reference matrix to identify nearest neighbors.

.. code-block:: bash

   beta_impute.py refknn \
       methylation.tsv \
       methylation.refknn.tsv \
       --reference reference.tsv \
       --neighbors 5

Random forest imputation
~~~~~~~~~~~~~~~~~~~~~~~~

Impute missing values using iterative random forest regression while
constraining predictions to the valid beta-value range.

.. code-block:: bash

   beta_impute.py rf \
       methylation.tsv \
       methylation.rf.tsv \
       --trees 500 \
       --max-iter 20 \
       --min-value 0 \
       --max-value 1

SoftImpute
~~~~~~~~~~

Perform low-rank matrix completion using SoftImpute.

.. code-block:: bash

   beta_impute.py softimpute \
       methylation.tsv \
       methylation.softimpute.tsv \
       --max-rank 20 \
       --max-iter 100 \
       --min-value 0 \
       --max-value 1

MOREL
~~~~~

Impute systematic block-wise missing values between two sample groups using
Random Forest.

.. code-block:: bash

   beta_impute.py morel \
       methylation.tsv \
       methylation.morel.tsv \
       --group groups.json \
       --model RF \
       --min-value 0 \
       --max-value 1

Genomic nearest-neighbor imputation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Impute missing values using nearby CpGs in the genome.

.. code-block:: bash

   beta_impute.py gnn \
       methylation.tsv \
       methylation.gnn.tsv \
       --gfile GRCh38_methyl_probes.info.tsv \
       --up-dist 200 \
       --down-dist 200 \
       --up-ncpg 3 \
       --down-ncpg 3

Evaluate imputation accuracy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare the imputed matrix with a complete reference ("truth") matrix and
report MAE, RMSE, and R² for the positions that were originally missing.

.. code-block:: bash

   beta_impute.py knn \
       methylation.missing.tsv \
       methylation.imputed.tsv \
       --neighbors 5 \
       --truth methylation.truth.tsv
