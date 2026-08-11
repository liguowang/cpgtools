beta_deconvolution
==================

Overview
--------

``beta_deconvolution`` estimates cell or tissue fractions from a DNA
methylation Beta-value matrix using a bundled reference methylation atlas.

The input matrix must have **CpGs in rows** and **samples in columns**. CpGs
shared between the input matrix and the selected reference atlas are used for
deconvolution.

Two reference atlases are available:

* ``CellTrace`` -- default reference atlas
* ``EpiDISH`` -- alternative reference atlas

For each sample, estimated component fractions are constrained to be
non-negative and are normalized to sum to 1.


Input
-----

The input is a tabular Beta-value matrix. The first column contains CpG IDs
and the remaining columns contain samples. Common delimiters are detected
automatically and compressed input is supported.

Example::

   CpG_ID      Sample_01   Sample_02   Sample_03
   cg00000029  0.432       0.517       0.481
   cg00000108  0.821       0.793       0.806
   cg00000165  0.135       0.164       0.142

Duplicate CpG IDs are reduced to the first occurrence. Non-numeric values are
treated as missing values, and only complete CpGs are used for fitting each
sample.


Cell Types
----------

The CellTrace reference panel includes immune-cell populations such as:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Label
     - Cell type
   * - ``Bmem``
     - B memory cells
   * - ``Bnv``
     - B naïve cells
   * - ``CD4Tmem``
     - CD4+ memory T cells
   * - ``CD4Tnv``
     - CD4+ naïve T cells
   * - ``CD8T``
     - CD8+ T cells
   * - ``CD8Tmem``
     - CD8+ memory T cells
   * - ``CD8Tnv``
     - CD8+ naïve T cells
   * - ``Treg``
     - Regulatory T cells
   * - ``NK``
     - Natural killer cells
   * - ``Bas``
     - Basophils
   * - ``Neu``
     - Neutrophils
   * - ``Eos``
     - Eosinophils
   * - ``Mono``
     - Monocytes

When corresponding subtype labels are present in the reference atlas,
``beta_deconvolution`` also creates a combined result table by collapsing:

* B memory + B naïve into ``Bcell.pred``
* CD4 memory + CD4 naïve into ``CD4T.pred``
* CD8 memory + CD8 naïve into ``CD8T.pred``
* neutrophils + eosinophils + basophils into ``Granulocyte.pred``


Fitting Methods
---------------

Two fitting methods are available:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Method
     - Description
   * - ``nnls``
     - Non-negative least squares. Coefficients are estimated with non-negative
       constraints and then normalized to sum to 1.
   * - ``bounded``
     - Bounded least squares with coefficients constrained to [0, 1], followed
       by normalization to sum to 1.


Usage
-----

Basic usage::

   beta_deconvolution samples.tsv

Select a reference atlas and fitting method::

   beta_deconvolution samples.tsv \
       --atlas CellTrace \
       --method nnls \
       -o output

Useful options include:

* ``--atlas {CellTrace,EpiDISH}`` -- reference atlas
* ``--method {nnls,bounded}`` -- fitting method
* ``-p``, ``--processes`` -- number of worker processes
* ``-r``, ``--residuals`` -- write per-sample fit metrics
* ``--plot`` -- generate a stacked-bar plot
* ``--slim`` -- write result values without row or column labels
* ``-o``, ``--out_prefix`` -- output prefix

If no output prefix is supplied, the input filename without its extension is
used.

Display all options with::

   beta_deconvolution -h


Output
------

For output prefix ``output``, the main result is:

* ``output_deconv.tsv`` -- estimated fractions for each sample

Additional files may include:

* ``output_combined_deconv.tsv`` -- immune subtypes collapsed into broader
  categories when applicable
* ``output_fit_metrics.tsv`` -- L2 residual, RMSE, and number of CpGs used for
  each sample; generated with ``--residuals``
* ``output_deconv_plot.png`` -- stacked-bar plot of estimated fractions;
  generated with ``--plot``

The main result table contains samples in rows and estimated cell/tissue
fractions in columns.


Fit Metrics
-----------

When ``--residuals`` is used, the following metrics are reported for each
sample:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Metric
     - Description
   * - ``residual_l2``
     - L2 norm of the residual vector between observed and reconstructed
       methylation values.
   * - ``rmse``
     - Root mean squared error of the fitted methylation profile.
   * - ``n_cpg``
     - Number of complete CpGs used for fitting the sample.


Plot
----

With ``--plot``, a stacked-bar plot of estimated fractions is generated.
Components with estimated fractions below 1% are grouped into an ``other``
category for visualization.
