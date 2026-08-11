CpG_anno_position
=================

Overview
--------

``CpG_anno_position`` annotates CpG sites or other BED intervals by overlap
with user-defined genomic regions.

Both the input CpG file and annotation file must contain at least three BED
columns: chromosome, start, and end. Additional columns in the CpG file are
preserved in the output.

If annotation column 4 is present, it is used as the annotation label.
Otherwise, a label is generated from the genomic coordinates, for example
``chr1:100-200``.

When multiple annotation regions overlap the same CpG, unique labels are
sorted and joined with a configurable separator.


Input Files
-----------

CpG file
~~~~~~~~

The input must be BED3 or BED3+. Compressed input is supported by the
CpGtools reader.

Example::

   chr1    10847    10848    cg26928153    0.8965    +
   chr1    10849    10850    cg16269199    0.7915    +
   chr1    15864    15865    cg13869341    0.9325    +

The CpG file and annotation file must use the same genome assembly.


Annotation file
~~~~~~~~~~~~~~~

The annotation file must be BED3 or BED3+.

Example::

   chr1    10000    12000    FOXA2
   chr1    15000    17000    MNT

If column 4 is absent, the annotation label is generated from the interval
coordinates.


Centered Annotation Window
--------------------------

By default, the full annotation interval is used.

With ``--window N``, each annotation region is restricted to an
``N``-base-pair window centered on the original interval midpoint. The window
is clipped to the original annotation interval and therefore never extends
beyond it.

For example::

   CpG_anno_position \
       -i cpg.bed \
       -a annotations.bed \
       --window 1000 \
       -o output

If ``--window 0`` is used, the full annotation interval is retained.


Usage
-----

Basic usage::

   CpG_anno_position \
       -i test_01.hg19.bed6 \
       -a hg19_ENCODE_338TF_130Cell_E3.bed.gz \
       -o output

If the first non-directive line of the CpG input is a header::

   CpG_anno_position \
       -i test_01.hg19.bed6 \
       -a hg19_ENCODE_338TF_130Cell_E3.bed.gz \
       --header \
       -o output

Useful options include:

* ``-w``, ``--window`` -- centered annotation-window size in bp
  (default: 0)
* ``--separator`` -- separator used when multiple labels overlap the same CpG
  (default: ``,``,)
* ``--na_rep`` -- value written when no annotation overlaps a CpG
  (default: ``NA``)
* ``--annotation_name`` -- name of the appended annotation column when
  ``--header`` is used; by default, the annotation filename is used
* ``-l``, ``--header`` -- treat the first non-directive CpG line as a header
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   CpG_anno_position -h


Output
------

For output prefix ``output``, the command writes:

* ``output.anno.tsv``

The original CpG columns are preserved and one annotation column is appended.

Example::

   #Chrom   Start   End     Name         Beta     Strand   Annotation
   chr1     10847   10848   cg26928153   0.8965   +        NA
   chr1     10849   10850   cg16269199   0.7915   +        NA
   chr1     15864   15865   cg13869341   0.9325   +        NA
   chr1     534241  534242  cg24669183   0.7941   +        FOXA2,MNT

CpGs without an overlapping annotation receive the value specified by
``--na_rep``.

Invalid CpG records are preserved in the output and receive the same missing
annotation value.


Pre-computed Annotation Datasets
--------------------------------

Several hg19 annotation datasets are available from the CpGtools project,
including:

* `ENCODE transcription-factor binding sites
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_338TF_130Cell_E3.bed.gz>`_
* `ENCODE DNase I hypersensitive sites
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_DNaseI_125Cells_V3.bed.gz>`_
* `ENCODE ChromHMM states
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_chromHMM_states_9Cells.merge.bed.gz>`_
* `FANTOM5 enhancers
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_FANTOM_enhancers_phase_1_and_2.bed.gz>`_
* `ENCODE H3K4me1 peaks
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_H3K4me1_11_cellLines_ChIP.bed.gz>`_
* `ENCODE H3K4me3 peaks
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_H3K4me3_11_cellLines_ChIP.bed.gz>`_
* `ENCODE H3K27ac peaks
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_H3K27ac_11_cellLines_ChIP.bed.gz>`_

hg38/GRCh38 versions are available from the
`CpGtools data directory
<https://sourceforge.net/projects/cpgtools/files/data/>`_.


Example Data
------------

* `CpG BED6 file
  <https://sourceforge.net/projects/cpgtools/files/test/test_01.hg19.bed6>`_
* `ENCODE TF annotation file
  <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_338TF_130Cell_E3.bed.gz/download>`_
