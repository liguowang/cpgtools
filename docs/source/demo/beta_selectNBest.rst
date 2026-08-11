beta_selectNBest
================

Overview
--------

``beta_selectNBest`` selects CpG features associated with a categorical
phenotype using scikit-learn feature-selection methods.

The input matrix must have **CpGs in rows** and **samples in columns**.
CpGs containing missing values are removed before feature selection, and
zero-variance CpGs are removed by default.


Scoring Methods
---------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Method
     - Description
   * - ``anova``
     - ANOVA F-statistic for association between methylation values and group
       labels.
   * - ``mi``
     - Mutual information between methylation values and group labels.
   * - ``chisq``
     - Chi-square statistic. Input values must be non-negative.


Input Files
-----------

Beta matrix
~~~~~~~~~~~

A tab-delimited matrix with CpG IDs in the first column and samples in the
remaining columns. Compressed input is supported.

Example::

   CpG_ID   Sample_01   Sample_02   Sample_03   Sample_04
   cg_001   0.831035    0.878022    0.794427    0.880911
   cg_002   0.249544    0.209949    0.234294    0.236680
   cg_003   0.845065    0.843957    0.840184    0.824286


Group file
~~~~~~~~~~

A two-column CSV or TSV file containing sample IDs and categorical group
labels. A header is optional.

Example::

   Sample,Group
   Sample_01,normal
   Sample_02,normal
   Sample_03,tumor
   Sample_04,tumor


Usage
-----

Basic usage::

   beta_selectNBest \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp.csv \
       -k 100 \
       -s chisq \
       -o selected

Useful options include:

* ``-k``, ``-c``, ``--top``, ``--topK`` -- number of CpGs to select
  (default: 100)
* ``-s``, ``--score_function``, ``--score-function`` -- ``anova``, ``mi``,
  or ``chisq`` (default: ``chisq``)
* ``--random_state`` -- random seed for mutual-information scoring
* ``--keep_constant`` -- retain zero-variance CpGs
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   beta_selectNBest -h


Output
------

For output prefix ``selected``, the command writes:

* ``selected.selected_features.tsv`` -- selected CpG-by-sample matrix
* ``selected.feature_scores.tsv`` -- scores, p-values, ranks, and selection
  status for all retained CpGs
* ``selected.selected_cpgs.txt`` -- selected CpG IDs


Example Data
------------

* `Beta matrix <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
