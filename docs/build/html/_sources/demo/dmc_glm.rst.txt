dmc_glm
=======

Overview
--------

``dmc_glm`` performs differential CpG analysis using a Gaussian generalized
linear model (GLM) fitted to DNA methylation Beta-values.

The input matrix must have **CpGs/probes in rows** and **samples in columns**.
The first variable in the group/covariate file is treated as the primary
grouping variable; additional variables are included as covariates.

For each CpG, the model fitted in R is::

   glm(y ~ covariates, family = gaussian)

Raw p-values for the primary grouping variable are extracted from the R model
output and adjusted for multiple testing using the CpGtools
Benjamini-Hochberg procedure.


Requirements
------------

``dmc_glm`` requires ``Rscript`` to be available in ``PATH``.


Input Files
-----------

Beta-value matrix
~~~~~~~~~~~~~~~~~

The first row contains sample IDs and the first column contains CpG/probe IDs.

Example::

   CpG_ID      Sample_01   Sample_02   Sample_03   Sample_04
   cg00001099  0.775       0.812       0.623       0.598
   cg00000363  0.611       0.602       0.470       0.455

Sample IDs must be unique.

Non-numeric values and ``NaN`` are treated as missing values. Plain-text,
``.gz``, and ``.bz2`` inputs supported by the CpGtools reader may be used.


Group and covariate file
~~~~~~~~~~~~~~~~~~~~~~~~

The group/covariate file is parsed by CpGtools and may contain a primary
grouping variable plus additional covariates.

The first variable must be categorical. Additional variables may be
categorical or continuous.

Example::

   Sample,Group,Sex,Age
   Sample_01,Control,F,52
   Sample_02,Control,M,47
   Sample_03,Case,F,61
   Sample_04,Case,M,58

Every sample listed in the group file must occur in the methylation matrix.

Samples present in the methylation matrix but absent from the group file are
ignored by assigning ``NA`` to their covariates before model fitting.


Usage
-----

Basic two-group analysis::

   dmc_glm \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp.csv \
       -o GLM_2G

Analysis with additional covariates::

   dmc_glm \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp2.csv \
       -o GLM_2G

Available options are:

* ``-i``, ``--input_file`` -- Beta-value matrix
* ``-g``, ``--group`` -- group/covariate file
* ``-o``, ``--output`` -- output prefix
* ``--version`` -- show the CpGtools version

Display all options with::

   dmc_glm -h


Output
------

For output prefix ``GLM_2G``, the command writes:

* ``GLM_2G.r`` -- generated R script;
* ``GLM_2G.results.txt`` -- full GLM coefficient and p-value table;
* ``GLM_2G.warnings.txt`` -- stderr produced by ``Rscript``;
* ``GLM_2G.pval.txt`` -- original methylation table with primary-variable
  p-values and adjusted p-values appended.


Final p-value table
~~~~~~~~~~~~~~~~~~~

``GLM_2G.pval.txt`` preserves the original input rows and appends:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Column
     - Description
   * - ``pval``
     - Raw p-value for the first coefficient associated with the primary
       grouping variable.
   * - ``adj.pval``
     - Benjamini-Hochberg FDR-adjusted p-value.

If a valid primary-variable p-value is unavailable for a CpG, both fields are
written as ``NaN``.


Full model results
~~~~~~~~~~~~~~~~~~

``GLM_2G.results.txt`` contains estimated coefficients and p-values for all
terms returned by the R GLM, including the intercept, primary grouping
variable, and any additional covariates.

For categorical variables, R creates factor-level coefficients. The exact
column names therefore depend on the variable names and factor levels.


Important Behavior
------------------

* CpGs with missing Beta-values are fitted using the available observations
  through R's ``na.omit`` behavior.
* Input rows with fewer values than the header are padded with ``NA``; extra
  values are ignored with a warning.
* The primary grouping variable must be categorical.
* When the primary grouping variable has multiple factor levels,
  ``dmc_glm`` currently extracts the **first** p-value column whose name starts
  with the primary variable name. The final ``pval`` therefore represents one
  factor coefficient rather than an omnibus test of the entire factor.
* If the per-CpG R model fails, that CpG may be absent from
  ``<prefix>.results.txt`` and receives ``NaN`` in the final p-value table.


Example Data
------------

* `Beta-value matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
* `Two-group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_
* `Group/covariate file
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp2.csv>`_
