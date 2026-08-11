dmc_fisher
==========

Overview
--------

``dmc_fisher`` performs two-group differential CpG analysis using
**Fisher's exact test** on methylated-read and total-read counts.

If a group contains multiple biological or technical replicates, counts are
summed across replicates before the test. Therefore, Fisher's exact test does
not model between-replicate variability.

The input matrix must have **CpGs in rows** and **samples in columns**.


Input Files
-----------

Methylation count matrix
~~~~~~~~~~~~~~~~~~~~~~~~

Each methylation value must be represented as::

   methylated_count,total_count

where both values are non-negative integers,
``methylated_count <= total_count``, and ``total_count > 0``.

Example::

   CpG_ID   Sample_01   Sample_02   Sample_03   Sample_04
   CpG_1    129,170     166,178     7,9          6,16
   CpG_2    24,77       67,99       0,85         1,37

The first row contains sample IDs and the first column contains CpG/probe IDs.
Sample IDs must be unique.

Plain text, compressed files, and other input types supported by the CpGtools
reader may be used.

Fields that do not match the ``methylated,total`` format are treated as
missing and ignored. Syntactically valid count pairs with zero total coverage
or methylated counts greater than total counts are treated as errors.


Group file
~~~~~~~~~~

The group file is a comma-separated, two-column file with a header:

* column 1: sample ID
* column 2: group ID

Exactly two groups are required.

Example::

   Sample,Group
   Sample_01,Control
   Sample_02,Control
   Sample_03,Case
   Sample_04,Case

Sample IDs in the group file must be unique and must occur in the methylation
matrix. Samples present in the methylation matrix but absent from the group
file are ignored.

The two group IDs are sorted internally before constructing the contingency
table.


Method
------

For each CpG, methylated and unmethylated reads are summed within each group.

The resulting 2 x 2 table is::

                     Methylated    Unmethylated
   Group 1               m1           n1 - m1
   Group 2               m2           n2 - m2

Fisher's exact test is then applied using ``scipy.stats.fisher_exact``.

P-values are adjusted for multiple testing using the CpGtools
``multiple_testing_correction`` implementation, which performs
Benjamini-Hochberg false-discovery-rate correction.


Usage
-----

Basic usage::

   dmc_fisher \
       -i test_09.tsv.gz \
       -g test_09.grp.csv \
       -o test_fisher

Available options are:

* ``-i``, ``--input_file`` -- methylated/total count matrix
* ``-g``, ``--group`` -- two-group sample file
* ``-o``, ``--output`` -- output prefix
* ``--version`` -- show the CpGtools version

Display all options with::

   dmc_fisher -h


Output
------

For output prefix ``test_fisher``, the command writes:

* ``test_fisher.pval.txt``

The original input table is preserved and three columns are appended:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Column
     - Description
   * - ``OddsRatio``
     - Fisher's exact-test odds ratio.
   * - ``pval``
     - Fisher's exact-test p-value.
   * - ``adj.pval``
     - Benjamini-Hochberg FDR-adjusted p-value.

Example::

   ID               Sample_01   Sample_02   Sample_03   Sample_04   OddsRatio   pval       adj.pval
   chr10:100011340  12,14       26,37       0,18        10,24       9.353846    1.21e-06   6.34e-05
   chr10:100011388  18,18       47,54       19,23       18,19       1.254826    0.7574     1.0

If Fisher's exact test returns a missing statistic, it is written as ``NA``.


Notes
-----

* Exactly two groups are required.
* Replicate counts are pooled within each group; replicate-level variability is
  not modeled.
* Input rows with fewer data fields than the header are padded with missing
  values; extra fields are ignored with a warning.
* Missing statistical results are excluded from multiple-testing correction
  and retained as missing in the output.


Example Data
------------

* `Count matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_09.tsv.gz/download>`_
* `Group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_09.grp.csv/download>`_
