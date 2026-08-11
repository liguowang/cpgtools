dmc_logit
=========

Overview
--------

``dmc_logit`` performs differential CpG analysis using logistic regression on
methylated-read and total-read counts.

The input matrix must have **CpGs/probes in rows** and **samples in columns**.
Each methylation value is represented as::

   methylated_count,total_count

where both values are non-negative integers,
``methylated_count <= total_count``, and ``total_count > 0``.

The first variable in the group/covariate file is treated as the primary
grouping variable. Additional variables are included as covariates.


Statistical Model
-----------------

For each CpG, the model fitted in R is::

   glm(cbind(m, t - m) ~ covariates, family = ...)

where ``m`` is the methylated-read count and ``t`` is the total-read count.

Two GLM families are available:

.. list-table::
   :header-rows: 1
   :widths: 18 28 54

   * - Value
     - Family
     - Description
   * - ``1``
     - ``quasibinomial``
     - Default. Allows for overdispersion relative to a standard binomial
       model.
   * - ``2``
     - ``binomial``
     - Standard binomial logistic regression.

Raw p-values for the primary grouping variable are extracted from the R model
output and adjusted using the CpGtools Benjamini-Hochberg multiple-testing
procedure.


Requirements
------------

``dmc_logit`` requires ``Rscript`` to be available in ``PATH``.


Input Files
-----------

Methylation count matrix
~~~~~~~~~~~~~~~~~~~~~~~~

The first row contains sample IDs and the first column contains CpG/probe IDs.

Example::

   CpG_ID   Sample_01   Sample_02   Sample_03   Sample_04
   cg_001   12,14       26,37       0,18        10,24
   cg_002   18,18       47,54       19,23       18,19

Sample IDs must be unique.

Plain text, ``.gz``, ``.bz2``, and other formats supported by the CpGtools
reader may be used.

Fields that do not match the ``methylated,total`` pattern are treated as
missing values. Syntactically valid count pairs are considered errors when
``total_count <= 0`` or ``methylated_count > total_count``.


Group and Covariate File
~~~~~~~~~~~~~~~~~~~~~~~~

The group/covariate file may contain a primary grouping variable plus
additional covariates.

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
ignored by assigning ``NA`` to their count and covariate values before model
fitting.


Usage
-----

Quasibinomial model::

   dmc_logit \
       -i test_04_TwoGroup.tsv.gz \
       -g test_04_TwoGroup.grp.csv \
       -o output_quasibin

Binomial model::

   dmc_logit \
       -i test_04_TwoGroup.tsv.gz \
       -g test_04_TwoGroup.grp.csv \
       -f 2 \
       -o output_bin

Available options are:

* ``-i``, ``--input_file`` -- methylated/total count matrix
* ``-g``, ``--group`` -- group/covariate file
* ``-f``, ``--family {1,2}`` -- GLM family; 1 = ``quasibinomial``,
  2 = ``binomial`` (default: 1)
* ``-o``, ``--output`` -- output prefix
* ``--version`` -- show the CpGtools version

Display all options with::

   dmc_logit -h


Output
------

For output prefix ``output_quasibin``, the command writes:

* ``output_quasibin.r`` -- generated R script;
* ``output_quasibin.results.txt`` -- full model coefficient and p-value table;
* ``output_quasibin.warnings.txt`` -- stderr produced by ``Rscript``;
* ``output_quasibin.pval.txt`` -- original methylation table with
  primary-variable p-values and adjusted p-values appended.


Final p-value Table
~~~~~~~~~~~~~~~~~~~

``<prefix>.pval.txt`` preserves the original input rows and appends:

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


Full Model Results
~~~~~~~~~~~~~~~~~~

``<prefix>.results.txt`` contains estimated coefficients and p-values for all
terms returned by the R GLM, including the intercept, primary grouping
variable, and additional covariates.

For categorical variables, R creates factor-level coefficients. The exact
column names therefore depend on the variable names and factor levels.


Important Behavior
------------------

* Per-CpG fitting uses ``na.omit`` in R, so samples with missing count or
  covariate values are removed for that CpG.
* Input rows with fewer data fields than the header are padded with ``NA``;
  extra values are ignored with a warning.
* The primary grouping variable must be categorical.
* If all methylated counts for a CpG sum to zero, the model p-values are set
  to missing.
* When the primary grouping variable has multiple factor levels,
  ``dmc_logit`` currently extracts the **first** p-value column whose name
  starts with the primary variable name. The final ``pval`` therefore
  represents one factor coefficient rather than an omnibus test of the entire
  factor.
* If the per-CpG R model fails, that CpG may be absent from
  ``<prefix>.results.txt`` and receives ``NaN`` in the final p-value table.


Example Data
------------

* `Count matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_04_TwoGroup.tsv.gz>`_
* `Group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_04_TwoGroup.grp.csv>`_
