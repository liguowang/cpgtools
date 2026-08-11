dmc_nonparametric
=================

Overview
--------

``dmc_nonparametric`` performs differential CpG analysis using nonparametric
tests on DNA methylation Beta-values.

The test is selected automatically from the number of biological groups:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Number of groups
     - Test
   * - 2
     - Two-sided Mann-Whitney U test
   * - 3 or more
     - Kruskal-Wallis H-test

The input matrix must have **CpGs/probes in rows** and **samples in columns**.


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

Non-numeric and non-finite values are treated as missing and ignored.
Plain text, ``.gz``, ``.bz2``, and other formats supported by the CpGtools
reader may be used.


Group file
~~~~~~~~~~

The group file is a comma-separated, two-column file with a header:

* column 1: sample ID
* column 2: group ID

At least two groups are required.

Example::

   Sample,Group
   Sample_01,Control
   Sample_02,Control
   Sample_03,Case
   Sample_04,Case

Sample IDs in the group file must be unique and must occur in the
methylation matrix.

Samples present in the methylation matrix but absent from the group file are
ignored.


Method
------

Two groups
~~~~~~~~~~

For exactly two groups, ``dmc_nonparametric`` performs a two-sided
Mann-Whitney U test using ``scipy.stats.mannwhitneyu``.

If either group has no valid observations for a CpG, the test result is
missing.


Three or more groups
~~~~~~~~~~~~~~~~~~~~

For three or more groups, the command performs a Kruskal-Wallis H-test using
``scipy.stats.kruskal``.

Groups with no valid observations for a particular CpG are omitted from that
CpG's test. If fewer than two groups contain valid observations, the result is
missing.

After testing all CpGs, valid p-values are adjusted using the CpGtools
Benjamini-Hochberg multiple-testing procedure.


Usage
-----

Two-group analysis::

   dmc_nonparametric \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp.csv \
       -o U_test

Three-group analysis::

   dmc_nonparametric \
       -i test_06_ThreeGroup.tsv.gz \
       -g test_06_ThreeGroup.grp.csv \
       -o H_test

Available options are:

* ``-i``, ``--input_file`` -- Beta-value matrix
* ``-g``, ``--group`` -- sample/group file
* ``-o``, ``--output`` -- output prefix
* ``--version`` -- show the CpGtools version

Display all options with::

   dmc_nonparametric -h


Output
------

For output prefix ``U_test``, the command writes:

* ``U_test.pval.txt``

The original input table is preserved and two columns are appended:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Column
     - Description
   * - ``pval``
     - Raw Mann-Whitney U or Kruskal-Wallis p-value.
   * - ``adj.pval``
     - Benjamini-Hochberg FDR-adjusted p-value.

Missing statistical results are written as ``NaN`` and are excluded from
multiple-testing correction.


Input-row Handling
------------------

If a CpG row contains fewer data fields than the header, missing values are
padded with ``NaN``. Extra fields are ignored. A warning is written in either
case.

The original input line is preserved in the final output; only ``pval`` and
``adj.pval`` are appended.


Example Data
------------

Two groups:

* `Beta-value matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
* `Group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_

Three groups:

* `Beta-value matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_06_ThreeGroup.tsv.gz>`_
* `Group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_06_ThreeGroup.grp.csv>`_
