beta_stats
==========

Overview
--------

``beta_stats`` summarizes CpG methylation values within user-defined genomic
regions.

For each region, six columns are appended:

* ``CpG_count``
* ``Minimum_beta``
* ``Maximum_beta``
* ``Mean_beta``
* ``Median_beta``
* ``Standard_deviation``

BED coordinates are interpreted as **0-based, half-open** intervals.


Input Files
-----------

CpG methylation file
~~~~~~~~~~~~~~~~~~~~

A BED6 or BED6+ file. Column 5 must contain Beta-values. Compressed input is
supported.

Example::

   chr22   44021512    44021513    cg24055475    0.9231    -
   chr13   111568382   111568383   cg06540715    0.1071    +
   chr20   44033594    44033595    cg21482942    0.6122    -


Region file
~~~~~~~~~~~

A BED3 or BED3+ file containing the genomic regions to summarize. All original
region columns are preserved in the output.

Example::

   chr1    15864    15865
   chr1    18826    18827
   chr1    29406    29407


Usage
-----

Basic usage::

   beta_stats \
       -i test_02.bed6.gz \
       -r hg19.RefSeq.union.1Kpromoter.bed.gz \
       -o region_stats

Useful options include:

* ``--header`` -- add a header row to the output
* ``--na_rep`` -- text used when statistics are unavailable (default: ``NA``)
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   beta_stats -h


Output
------

For output prefix ``region_stats``, the command writes:

* ``region_stats.region_stats.tsv``

Invalid region records are retained and filled with the value specified by
``--na_rep``.


Example Data
------------

* `CpG BED6 file <https://sourceforge.net/projects/cpgtools/files/test/test_02.bed6.gz>`_
* `Example promoter regions <https://sourceforge.net/projects/cpgtools/files/test/hg19.RefSeq.union.1Kpromoter.bed.gz/download>`_
