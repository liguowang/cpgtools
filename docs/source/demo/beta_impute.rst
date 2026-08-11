beta_impute
===========

Overview
--------

``beta_impute`` provides a unified command-line interface for inspecting,
simulating, and imputing missing values in DNA methylation matrices.

The input matrix is expected to have **CpGs in rows** and **samples in
columns**. Common delimiters and compressed input files are detected
automatically.

Most imputation commands can optionally evaluate imputed values against a
truth matrix using MAE, RMSE, and :math:`R^2`.


Commands
--------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Command
     - Description
   * - ``toy``
     - Generate a synthetic methylation matrix with missing values.
   * - ``insertna``
     - Insert missing values until a requested total is reached.
   * - ``countna``
     - Report missing-value counts by CpG and sample.
   * - ``dropna``
     - Remove rows or columns containing missing values.
   * - ``constant``
     - Replace missing values with a fixed value.
   * - ``mean``
     - Impute using row- or column-wise means.
   * - ``median``
     - Impute using row- or column-wise medians.
   * - ``min`` / ``max``
     - Impute using row- or column-wise minima or maxima.
   * - ``rand``
     - Impute using randomly selected observed values from the same row or
       column.
   * - ``refknn``
     - KNN imputation using an external complete reference matrix.
   * - ``mw``
     - Impute from neighboring values in a moving window.
   * - ``knn``
     - KNN imputation within the input matrix.
   * - ``buck``
     - Iterative linear-regression imputation based on Buck's method.
   * - ``rf``
     - Iterative Random Forest regression.
   * - ``softimpute``
     - Low-rank matrix completion using SoftImpute.
   * - ``morel``
     - Impute systematic block-wise missingness using Random Forest or a dense
       neural network, with KNN for sporadic missingness.
   * - ``gnn``
     - Impute using genomically neighboring CpGs.


General Usage
-------------

Display all subcommands::

   beta_impute -h

Display options for a specific method::

   beta_impute knn -h

Most imputation commands use two positional arguments::

   beta_impute METHOD input_file output_file

Common options include:

* ``--decimal`` -- number of decimal places written to output (default: 5)
* ``--overwrite`` -- replace an existing output file
* ``--truth FILE`` -- evaluate imputation accuracy where supported

Enable debug logging at the top level with::

   beta_impute --debug knn input.tsv output.tsv


Examples
--------

KNN imputation
~~~~~~~~~~~~~~

By default, KNN searches for neighboring samples, uses three neighbors, and
weights neighbors by distance.

::

   beta_impute knn \
       methylation.tsv \
       methylation.knn.tsv \
       --neighbors 5 \
       --weights distance

Use ``--axis index`` to search neighboring CpGs instead of samples.


Reference-based KNN
~~~~~~~~~~~~~~~~~~~

The external reference matrix must not contain missing values.

::

   beta_impute refknn \
       methylation.tsv \
       methylation.refknn.tsv \
       --reference reference.tsv \
       --search_axis columns \
       --neighbors 3


Random Forest
~~~~~~~~~~~~~

::

   beta_impute rf \
       methylation.tsv \
       methylation.rf.tsv \
       --trees 500 \
       --max-iter 20

By default, Random Forest predictions are constrained to the interval
[0, 1] using ``--min-value 0`` and ``--max-value 1``.


SoftImpute
~~~~~~~~~~

::

   beta_impute softimpute \
       methylation.tsv \
       methylation.softimpute.tsv \
       --max-rank 20 \
       --max-iter 100

SoftImpute does not apply lower or upper bounds by default. For Beta-value
matrices, add ``--min-value 0`` and ``--max-value 1`` if bounded output is
desired.


Evaluating Imputation
---------------------

Most imputation commands accept ``--truth FILE``.

Evaluation is restricted to cells that:

* were missing in the original input;
* are observed in the truth matrix; and
* are present after imputation.

The truth matrix must contain all CpG IDs and sample IDs from the original
input.

Reported metrics are:

* **MAE** -- mean absolute error
* **RMSE** -- root mean squared error
* **R2** -- coefficient of determination

Example::

   beta_impute knn \
       methylation.tsv \
       methylation.knn.tsv \
       --truth truth.tsv


Missing-value Utilities
-----------------------

Generate a toy matrix
~~~~~~~~~~~~~~~~~~~~~

::

   beta_impute toy toy.tsv \
       --rows 1000 \
       --cols 20 \
       --missingness 0.1

``--missingness`` values between 0 and 1 specify a fraction; values greater
than or equal to 1 specify an approximate number of missing entries.


Insert missing values
~~~~~~~~~~~~~~~~~~~~~

``target_missing`` is the desired **total** number of missing values in the
output, including any missing values already present.

::

   beta_impute insertna \
       methylation.tsv \
       methylation.missing.tsv \
       5000


Count missing values
~~~~~~~~~~~~~~~~~~~~

::

   beta_impute countna methylation.tsv

By default, the command writes:

* ``missing_by_cpg.tsv``
* ``missing_by_sample.tsv``

Use ``--row-report`` and ``--column-report`` to change these filenames.


Drop incomplete rows or columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, incomplete CpG rows are removed.

::

   beta_impute dropna \
       methylation.tsv \
       methylation.complete.tsv

To remove incomplete samples instead::

   beta_impute dropna \
       methylation.tsv \
       methylation.complete.tsv \
       --axis columns


MOREL
-----

``morel`` targets systematic block-wise missingness between two sample groups.
Sporadic missing values are first handled by KNN, and block-wise missing values
are then predicted using either Random Forest (``RF``) or a dense neural
network (``DNN``).

If ``--group`` is omitted, the two sample groups are inferred from missingness
patterns using K-means clustering.

When ``--group`` is supplied, the current implementation expects a **JSON**
object mapping two group names to sample-name lists, for example::

   {
       "group_1": ["Sample_01", "Sample_02"],
       "group_2": ["Sample_03", "Sample_04"]
   }

Example::

   beta_impute morel \
       methylation.tsv \
       methylation.morel.tsv \
       --group groups.json \
       --model RF

Important defaults include:

* ``--knn-neighbors 5``
* ``--knn-weights uniform``
* ``--n-iter 10``
* ``--n-estimators 100``
* ``--max-depth 30``
* ``--min-value 0``
* ``--max-value 1``

TensorFlow is required only when ``--model DNN`` is used.


Genomic Nearest-neighbor Imputation
-----------------------------------

``gnn`` estimates missing values using nearby CpGs. A genomic annotation file
is required and must contain at least five columns:

``chrom``, ``start``, ``end``, ``cpg_id``, and ``CRE``.

Example::

   beta_impute gnn \
       methylation.tsv \
       methylation.gnn.tsv \
       --gfile GRCh38_methyl_probes.info.tsv \
       --up-dist 200 \
       --down-dist 200

Important options include:

* ``--up-dist`` / ``--down-dist`` -- maximum search distance in bp
  (default: 100)
* ``--up-ncpg`` / ``--down-ncpg`` -- maximum neighboring CpGs used on each
  side (default: 2)
* ``--same-cre`` -- require neighbors to share a candidate cis-regulatory
  element with the target CpG
* ``--method {AA,WA,TA}`` -- arithmetic average, inverse-distance weighted
  average, or trimmed average (default: ``TA``)
* ``--cpgfile`` -- append additional CpG IDs to the matrix for attempted
  imputation

Display all GNN options with::

   beta_impute gnn -h
