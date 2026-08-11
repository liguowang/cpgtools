dmc_bb
======

Overview
--------

``dmc_bb`` performs differential CpG analysis with a **beta-binomial
regression model** using methylated-read and total-read counts.

The input matrix must have **CpGs in rows** and **samples in columns**. Each
methylation value is represented as ``c,n``, where:

* ``c`` is the number of methylated reads;
* ``n`` is the total number of reads;
* ``0 <= c <= n``; and
* ``n > 0``.

The model is fitted in R with the ``aod::betabin`` function and can include
categorical or continuous covariates.


Requirements
------------

``dmc_bb`` requires:

* ``Rscript`` available in ``PATH``;
* the R package `aod
  <https://cran.r-project.org/web/packages/aod/index.html>`_.

Install ``aod`` from R with::

   install.packages("aod")


Input Files
-----------

Methylation count matrix
~~~~~~~~~~~~~~~~~~~~~~~~

The first row contains sample IDs and the first column contains CpG/probe IDs.
Plain-text, ``.gz``, and ``.bz2`` input files are supported.

Example::

   CpG_ID   A_1      A_2      A_3      B_1      B_2      B_3
   CpG_1    129,170  166,178  7,9      6,16     10,10    10,15
   CpG_2    0,77     0,99     0,85     1,37     3,37     0,42

Values that do not match the ``c,n`` pattern are treated as missing values.
Valid count pairs must have positive total coverage and cannot have more
methylated reads than total reads.

This command expects **count/proportion data, not Beta-values**.


Group and covariate file
~~~~~~~~~~~~~~~~~~~~~~~~

The group file defines the biological group and optional covariates for each
sample.

The first variable is the grouping variable and must be categorical.
Additional variables may be categorical or continuous.

For example, a group file might contain::

   Sample,Group,Sex,Age
   A_1,A,F,45
   A_2,A,M,51
   A_3,A,F,48
   B_1,B,M,47
   B_2,B,F,52
   B_3,B,M,50

Sample IDs must correspond to those in the methylation matrix.


Model
-----

For each CpG, ``dmc_bb`` fits a beta-binomial regression of the form::

   cbind(methylated, total - methylated) ~ Group + covariates

with a logit link and a constant dispersion model.

Rows containing missing methylation values or missing covariates are omitted
by R for that CpG before model fitting.


Usage
-----

Basic usage::

   dmc_bb \
       -i test_04_TwoGroup.tsv.gz \
       -g test_04_TwoGroup.grp.csv \
       -o OUT_bb

Available options are:

* ``-i``, ``--input_file`` -- methylated/total count matrix
* ``-g``, ``--group`` -- group and covariate file
* ``-o``, ``--output`` -- output prefix
* ``--version`` -- show the CpGtools version

Display help with::

   dmc_bb -h


Output
------

For output prefix ``OUT_bb``, the command creates:

* ``OUT_bb.results.txt`` -- beta-binomial regression coefficients and p-values
  for each CpG;
* ``OUT_bb.r`` -- generated R script used for the analysis;
* ``OUT_bb.warnings.txt`` -- messages written by ``Rscript`` to standard error.

The result columns are generated from the fitted model terms. For each model
coefficient, the output contains its estimated coefficient and corresponding
Wald z-test p-value.

The first column is always::

   ID

and the remaining columns follow the model variables, for example::

   ID   (Intercept).coef   Group.coef   Age.coef   (Intercept).pval   Group.pval   Age.pval

The exact coefficient names depend on the grouping variable, covariates, and
their factor levels.


Notes
-----

* ``dmc_bb`` supports more than two biological groups when the grouping
  variable contains multiple factor levels.
* Covariates are included directly in the beta-binomial regression model.
* CpGs with no methylated reads across the samples used by the model receive
  missing p-values.
* Because the analysis is performed through a generated R script, inspect
  ``<prefix>.warnings.txt`` if a model fails or R reports warnings.


Example Data
------------

* `Count matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_04_TwoGroup.tsv.gz>`_
* `Group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_04_TwoGroup.grp.csv>`_
