CpG_aggregation
===============

Overview
--------

``CpG_aggregation`` aggregates CpG methylation values over user-defined
genomic regions such as promoters, CpG islands, exons, or other BED intervals.

Two input data types are supported:

* ``count`` -- methylated and total read counts, for example ``3,10``
* ``beta`` -- Beta-values, for example ``0.30``

The CpG methylation value must be stored in **column 4** of the CpG BED file.


Input Files
-----------

CpG file
~~~~~~~~

The CpG input must be BED3+ with the methylation value in column 4.
Compressed input is supported.

Count-data example::

   chr1    100017748    100017749    3,10
   chr1    100017769    100017770    0,10
   chr1    100017853    100017854    16,21

For ``count`` input, values must satisfy:

* methylated reads >= 0
* total reads >= 0
* methylated reads <= total reads

Beta-value example::

   chr1    100017748    100017749    0.30
   chr1    100017769    100017770    0.05
   chr1    100017853    100017854    0.76

For ``beta`` input, finite values must fall within [0, 1].


Region file
~~~~~~~~~~~

The region file must be BED3+ and may contain any number of additional
columns. All original region columns are preserved in the output.

Example::

   chr1    567292    568293
   chr1    713567    714568
   chr1    762401    763402


Count-mode Outlier Filtering
----------------------------

In ``count`` mode, CpG-level outliers are filtered by default.

For each region, the original regional methylation proportion is calculated
from all overlapping CpGs with total read count greater than zero:

.. math::

   \beta = \frac{\sum m_i}{\sum n_i}

where :math:`m_i` is the methylated-read count and :math:`n_i` is the total
read count for CpG :math:`i`.

Each CpG is then tested against this regional methylation probability using a
**two-sided exact binomial test**. A Bonferroni-adjusted cutoff is used:

.. math::

   p_{\mathrm{cutoff}} = \frac{\alpha}{N_{\mathrm{CpG}}}

CpGs with :math:`p < p_{\mathrm{cutoff}}` are removed before the filtered
regional methylation value is calculated.

Use ``--no_outlier_filter`` to disable this filtering.

Outlier filtering is applied only in ``count`` mode.


Usage
-----

Count data::

   CpG_aggregation \
       -i test_03_RRBS.bed.gz \
       -b hg19.RefSeq.union.1Kpromoter.bed.gz \
       -t count \
       -o out

Beta-values::

   CpG_aggregation \
       -i methylation.beta.bed.gz \
       -b regions.bed.gz \
       -t beta \
       -o out

Useful options include:

* ``-t``, ``--type {count,beta}`` -- methylation data type in column 4;
  required
* ``-a``, ``--alpha`` -- family-wise error rate for count-mode outlier
  filtering (default: 0.05)
* ``--no_outlier_filter`` -- disable count-mode outlier filtering
* ``--min_cpg`` -- minimum number of valid overlapping CpGs required for a
  region (default: 1)
* ``--ddof {0,1}`` -- degrees-of-freedom adjustment used for the standard
  deviation in Beta-value mode (default: 0)
* ``--header`` -- write a header row
* ``--na_rep`` -- text used when summary values are unavailable
  (default: ``NA``)
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   CpG_aggregation -h


Output
------

For output prefix ``out``, the command writes:

* ``out.aggregation.tsv``

All columns from the region file are retained, and summary columns are
appended.


Count-mode Output
~~~~~~~~~~~~~~~~~

The following columns are appended:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``N_CpG_filtered``
     - Number of CpGs retained after outlier filtering.
   * - ``N_methyl_filtered``
     - Total methylated reads after filtering.
   * - ``N_total_filtered``
     - Total reads after filtering.
   * - ``Beta_filtered``
     - Aggregated methylation proportion after filtering.
   * - ``N_CpG_original``
     - Number of valid overlapping CpGs before filtering.
   * - ``N_methyl_original``
     - Total methylated reads before filtering.
   * - ``N_total_original``
     - Total reads before filtering.
   * - ``Beta_original``
     - Aggregated methylation proportion before filtering.

CpGs with total read count equal to zero are not included in count-mode
aggregation.


Beta-mode Output
~~~~~~~~~~~~~~~~

The following columns are appended:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Description
   * - ``N_CpG``
     - Number of overlapping Beta-values used.
   * - ``Mean_beta``
     - Mean Beta-value.
   * - ``Median_beta``
     - Median Beta-value.
   * - ``Min_beta``
     - Minimum Beta-value.
   * - ``Max_beta``
     - Maximum Beta-value.
   * - ``Std_beta``
     - Standard deviation of Beta-values.


Header and Missing Values
-------------------------

By default, no header is written. With ``--header``, the original region
columns are named:

``chrom``, ``start``, ``end``, ``column_4``, ``column_5``, ...

followed by the aggregation columns.

If a region has fewer than ``--min_cpg`` valid overlapping CpGs, the appended
summary fields are written using ``--na_rep``.

Invalid region records or records with an inconsistent number of columns are
preserved in the output, with missing summary values appended.


Example Data
------------

* `RRBS count data <https://sourceforge.net/projects/cpgtools/files/test/test_03_RRBS.bed.gz>`_
* `Example promoter regions <https://sourceforge.net/projects/cpgtools/files/test/hg19.RefSeq.union.1Kpromoter.bed.gz>`_
