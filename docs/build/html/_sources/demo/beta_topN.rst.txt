beta_topN
=========

Overview
--------

``beta_topN`` ranks CpGs by a row-wise summary statistic and retains the top
N features for downstream analyses such as PCA or clustering.

The input matrix should have **CpGs in rows** and **samples in columns**.


Ranking Scores
--------------

Supported row-wise scores are:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Score
     - Description
   * - ``std``
     - Standard deviation across samples. This is the default.
   * - ``var``
     - Variance across samples.
   * - ``mean``
     - Mean across samples.
   * - ``median``
     - Median across samples.
   * - ``mad``
     - Median absolute deviation across samples.


Input
-----

The input is a tabular CpG-by-sample matrix. Delimiters are detected
automatically and compressed input is supported.

Example::

   CpG_ID   Sample_01   Sample_02   Sample_03   Sample_04
   cg_001   0.831035    0.878022    0.794427    0.880911
   cg_002   0.249544    0.209949    0.234294    0.236680
   cg_003   0.845065    0.843957    0.840184    0.824286


Usage
-----

Basic usage::

   beta_topN \
       -i test_05_TwoGroup.tsv.gz \
       -c 500 \
       -o test_05_TwoGroup

Useful options include:

* ``-c``, ``--count``, ``--top`` -- number of CpGs to retain
  (default: 1000)
* ``-s``, ``--score`` -- ``std``, ``var``, ``mean``, ``median``, or ``mad``
* ``--ascending`` -- rank from smallest to largest score
* ``--na_policy {drop,omit}`` -- remove incomplete CpGs or calculate from
  available values
* ``--min_valid`` -- minimum observed values required when
  ``--na_policy omit`` is used (default: 2)
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   beta_topN -h


Output
------

For output prefix ``test_05_TwoGroup``, the command writes:

* ``test_05_TwoGroup.ranked.tsv`` -- all retained CpGs with ``Score`` and
  ``Rank`` columns
* ``test_05_TwoGroup.topN.tsv`` -- top-ranked CpG-by-sample matrix
* ``test_05_TwoGroup.topN_ids.txt`` -- selected CpG IDs


Example Data
------------

* `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
