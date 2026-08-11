beta_m_conversion
=================

Overview
--------

``beta_m_conversion`` converts DNA methylation **Beta-values to M-values** or
**M-values to Beta-values**.

The input matrix should have **CpGs in rows** and **samples in columns**. The
first column is preserved as the CpG identifier column.


Conversion
----------

Beta-values and M-values are related by:

.. math::

   M = \log_{2}\left(\frac{\beta}{1-\beta}\right)

and:

.. math::

   \beta = \frac{2^M}{2^M + 1}

When converting Beta-values to M-values, Beta-values are clamped away from 0
and 1 using a small epsilon to avoid infinite M-values.


Input
-----

The input may be plain text or compressed. Common delimiters are detected
automatically. Supported compression includes formats such as ``.gz``,
``.bz2``, and ``.xz``.

Example Beta-value matrix::

   CpG_ID   Sample_01   Sample_02   Sample_03   Sample_04
   cg_001   0.831035    0.878022    0.794427    0.880911
   cg_002   0.249544    0.209949    0.234294    0.236680
   cg_003   0.845065    0.843957    0.840184    0.824286

Non-numeric values in sample columns are treated as missing values.

Example data:

* `test_08.tsv.gz
  <https://sourceforge.net/projects/cpgtools/files/test/test_08.tsv.gz/download>`_


Usage
-----

Convert Beta-values to M-values::

   beta_m_conversion \
       -i test_08.tsv.gz \
       -d Beta \
       -o test_08

Convert M-values to Beta-values::

   beta_m_conversion \
       -i test_08.m.tsv \
       -d M \
       -o test_08

Useful options include:

* ``-d``, ``--dtype {Beta,M,beta,m}`` -- input data type
  (default: ``Beta``)
* ``-e``, ``--epsilon`` -- clamp Beta-values to
  ``[epsilon, 1 - epsilon]`` before conversion to M-values
  (default: ``1e-6``)
* ``-o``, ``--out_prefix`` -- output filename prefix

Display all options with::

   beta_m_conversion -h


Output
------

The output is tab-delimited and preserves the input CpG-by-sample layout.

The output filename is determined by the input data type:

* Beta input -> ``<prefix>.m.tsv``
* M input -> ``<prefix>.beta.tsv``

For example, with ``-o test_08`` and ``-d Beta``, the output is::

   test_08.m.tsv

Numeric values are written with six decimal places.
