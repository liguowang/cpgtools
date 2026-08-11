CpG_anno_probe
==============

Overview
--------

``CpG_anno_probe`` appends probe-level annotations to a tabular input file by
matching Illumina DNA methylation probe IDs, such as ``cg00000029``.

The input file may contain arbitrary columns. One column must contain probe
IDs. The annotation file must contain a probe-ID column plus one or more
annotation columns.

All original input columns are preserved, and annotation columns are appended.


Input Files
-----------

Input table
~~~~~~~~~~~

The input may be a regular or compressed tabular file. Tabs are preferred;
whitespace-delimited files are also accepted.

Example::

   Chrom   Start   End   Probe_ID     Beta     Strand
   chr1    10847   10848 cg26928153   0.8965   +
   chr1    10849   10850 cg16269199   0.7915   +
   chr1    15864   15865 cg13869341   0.9325   +

Use ``-p`` / ``--probe_column`` to specify the **zero-based** column containing
probe IDs. The default is column 0.

If the input contains a header, use ``--header`` (or ``-l``).


Annotation table
~~~~~~~~~~~~~~~~

The annotation file must contain a probe-ID column and at least one annotation
column.

Example::

   Probe_ID      hg19_pos        hg38_pos        geneSymbol
   cg26928153    chr1:10847      chr1:14465      DDX11L1
   cg16269199    chr1:10849      chr1:14467      DDX11L1

By default, the probe ID is expected in column 0. Use
``--annotation_probe_column`` to select another zero-based column.

Annotation-file headers can be controlled with::

   --annotation_header auto
   --annotation_header yes
   --annotation_header no

The default, ``auto``, recognizes common probe-ID header names including
``probeID``, ``probe_id``, ``IlmnID``, ``cpg``, ``cpg_id``, and ``name``.

If the annotation file has no header, generic annotation-column names are
generated.


Duplicate Probe IDs
-------------------

Duplicate probe IDs in the annotation file are controlled by
``--duplicate_policy``:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Policy
     - Behavior
   * - ``error``
     - Stop when a duplicate probe ID is encountered. This is the default.
   * - ``first``
     - Keep the first annotation record.
   * - ``last``
     - Keep the last annotation record.


Usage
-----

Basic usage::

   CpG_anno_probe \
       -i test_01.hg19.bed6 \
       -a MethylationEPIC_CpGtools.tsv.gz \
       -p 3 \
       --header \
       -o output

Useful options include:

* ``-p``, ``--probe_column`` -- zero-based probe-ID column in the input
  (default: 0)
* ``--annotation_probe_column`` -- zero-based probe-ID column in the annotation
  file (default: 0)
* ``-l``, ``--header`` -- treat the first input line as a header
* ``--annotation_header {auto,yes,no}`` -- annotation-header handling
  (default: ``auto``)
* ``--na_rep`` -- value appended when a probe has no matching annotation
  (default: ``NA``)
* ``--duplicate_policy {first,last,error}`` -- duplicate annotation handling
  (default: ``error``)
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   CpG_anno_probe -h


Output
------

For output prefix ``output``, the command writes:

* ``output.anno.tsv``

The output contains all original input columns followed by the annotation
columns.

When an input probe is absent from the annotation table, each appended
annotation field is filled with ``--na_rep``.

If ``--header`` is used, annotation-column names from the annotation file are
appended to the input header. For headerless annotation files, generic names
such as ``annotation_1`` and ``annotation_2`` are used.

Input rows with an inconsistent number of columns are preserved and receive
missing annotation values.


Example Data
------------

* `Example CpG file
  <https://sourceforge.net/projects/cpgtools/files/test/test_01.hg19.bed6>`_
* `MethylationEPIC annotation table
  <https://sourceforge.net/projects/cpgtools/files/data/MethylationEPIC_CpGtools.tsv.gz>`_


Annotation Content
------------------

The exact annotation fields are determined by the annotation file supplied to
``CpG_anno_probe``. For example, the CpGtools MethylationEPIC annotation table
may include genomic positions, gene annotations, CpG-island annotations, SNP
information, cross-reactivity flags, and regulatory annotations.

``CpG_anno_probe`` does not require a fixed set or fixed number of annotation
columns; it appends whatever annotation fields are present in the selected
annotation table.
