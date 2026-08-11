dmc_ttest
=========

Overview
--------

``dmc_ttest`` performs differential CpG analysis on DNA methylation
Beta-values using parametric tests.

The test is selected from the number of groups and the requested two-group
test mode:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Comparison
     - Test
   * - Two groups
     - Student's two-sample t-test by default.
   * - Two groups with ``--welch``
     - Welch's two-sample t-test.
   * - Two groups with ``--paired``
     - Paired t-test.
   * - Three or more groups
     - One-way ANOVA.

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

Non-numeric and non-finite values are treated as missing values. Plain-text,
``.gz``, ``.bz2``, and other formats supported by the CpGtools reader may be
used.


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

Sample IDs in the group file must be unique.

Samples listed in the group file but absent from the methylation matrix are
reported and excluded from the analysis. Samples present in the matrix but
not listed in the group file are ignored.


Two-group Tests
---------------

Student's t-test
~~~~~~~~~~~~~~~~

With two groups and no additional test flag, ``dmc_ttest`` uses the standard
independent two-sample t-test with equal variances assumed.

For a CpG to be tested, each group must contain at least two valid
observations.


Welch's t-test
~~~~~~~~~~~~~~

Use ``--welch`` for an independent two-sample t-test that does not assume
equal variances::

   dmc_ttest \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp.csv \
       --welch \
       -o ttest_welch


Paired t-test
~~~~~~~~~~~~~

Use ``--paired`` when the two groups contain matched samples::

   dmc_ttest \
       -i paired_beta.tsv \
       -g paired_groups.csv \
       --paired \
       -o ttest_paired

Important behavior:

* exactly two groups are required;
* the groups must contain the same number of samples after matching to the
  methylation matrix;
* samples are paired by their **order within each group**;
* missing values are removed pairwise for each CpG; and
* at least two complete pairs are required for a valid test.

If both ``--paired`` and ``--welch`` are supplied, ``--welch`` is ignored.


Multi-group ANOVA
-----------------

When three or more groups are defined, ``dmc_ttest`` automatically performs a
one-way ANOVA using ``scipy.stats.f_oneway``.

Each group must contain at least one valid observation for a CpG. If any group
has no valid value for that CpG, the ANOVA result is reported as missing.


Delta Beta
----------

For two-group comparisons, ``delta_beta`` is calculated as:

.. math::

   \Delta\beta
   =
   \mathrm{mean}(\mathrm{group}_1)
   -
   \mathrm{mean}(\mathrm{group}_2)

The group IDs are sorted internally, so group 1 is the first group in sorted
order and group 2 is the second.

For three or more groups, ``delta_beta`` is not defined and is written as
``n/a``.


Usage
-----

Standard two-group t-test::

   dmc_ttest \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp.csv \
       -o ttest_2G

Three-group ANOVA::

   dmc_ttest \
       -i test_06_ThreeGroup.tsv.gz \
       -g test_06_ThreeGroup.grp.csv \
       -o ttest_3G

Useful options include:

* ``-p``, ``--paired`` -- use a paired t-test for exactly two groups
* ``-w``, ``--welch`` -- use Welch's t-test for exactly two independent groups
* ``-i``, ``--input_file`` -- Beta-value matrix
* ``-g``, ``--group`` -- sample/group file
* ``-o``, ``--output`` -- output prefix
* ``--version`` -- show the CpGtools version

Display all options with::

   dmc_ttest -h


Output
------

For output prefix ``ttest_2G``, the command writes:

* ``ttest_2G.pval.txt``

The original input table is preserved and three columns are appended:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Column
     - Description
   * - ``delta_beta``
     - Difference between the two group means. Written as ``n/a`` for
       multi-group ANOVA or when it cannot be calculated.
   * - ``pval``
     - Raw p-value from the selected t-test or one-way ANOVA.
   * - ``adj.pval``
     - Benjamini-Hochberg FDR-adjusted p-value.

Missing p-values and adjusted p-values are written as ``n/a`` and excluded
from multiple-testing correction.


Input-row Handling
------------------

If a CpG row contains fewer data fields than the header, missing values are
padded with ``NaN``. Extra fields are ignored. A warning is reported in either
case.

The original input row is preserved in the final output; only
``delta_beta``, ``pval``, and ``adj.pval`` are appended.


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
