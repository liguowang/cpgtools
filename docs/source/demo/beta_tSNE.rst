beta_tSNE
=========

Overview
--------

``beta_tSNE`` performs t-distributed stochastic neighbor embedding (t-SNE) on
DNA methylation samples.

The input matrix must have **CpGs in rows** and **samples in columns**. CpGs
are treated as features and samples as observations.

A sample-group file is optional and can be used to color samples in the plot.


Input Files
-----------

Beta matrix
~~~~~~~~~~~

The first column contains CpG IDs and the remaining columns contain sample
Beta-values. Common delimiters and compressed input are supported.

Example::

   CpG_ID   Sample_01   Sample_02   Sample_03   Sample_04
   cg_001   0.831035    0.878022    0.794427    0.880911
   cg_002   0.249544    0.209949    0.234294    0.236680
   cg_003   0.845065    0.843957    0.840184    0.824286

CpG IDs and sample IDs must be unique.


Group file
~~~~~~~~~~

The optional group file contains two columns: sample ID and group label.
Comma- and tab-delimited files are supported, with or without a header.

Example::

   Sample,Group
   Sample_01,normal
   Sample_02,normal
   Sample_03,tumor
   Sample_04,tumor

If a group file is supplied, every sample in the Beta matrix must be present
in the group file.

If no group file is supplied, all samples are assigned to one group named
``All_samples``.


Missing Values and Preprocessing
--------------------------------

By default:

* CpGs containing any missing value are removed.
* zero-variance CpGs are removed.
* CpG features are standardized across samples before t-SNE.

Use ``--na_policy impute`` to replace missing values with the per-CpG mean.
CpGs containing only missing values are removed before imputation.

Use ``--no_standardize`` to skip feature standardization.


Usage
-----

Basic usage::

   beta_tSNE \
       -i cirrHCV_vs_normal.data.tsv \
       -g cirrHCV_vs_normal.grp.csv \
       -o HCV_vs_normal

Useful options include:

* ``-p``, ``--perplexity`` -- t-SNE perplexity (default: 5.0); must be
  smaller than the number of samples
* ``-n``, ``--n_components``, ``--ncomponent`` -- embedding dimensions,
  either 2 or 3 (default: 2)
* ``--max_iter``, ``--n_iter`` -- maximum optimization iterations
  (default: 5000)
* ``--learning_rate`` -- ``auto`` or a positive numeric value
  (default: ``auto``)
* ``--init {pca,random}`` -- initialization method (default: ``pca``)
* ``--metric`` -- distance metric accepted by scikit-learn t-SNE
  (default: ``euclidean``)
* ``--early_exaggeration`` -- early-exaggeration factor (default: 12.0)
* ``--seed`` -- random seed (default: 0)
* ``--na_policy {drop,impute}`` -- missing-value handling
  (default: ``drop``)
* ``--no_standardize`` -- skip CpG standardization
* ``--label`` -- add sample IDs to the plot
* ``--marker {circle,dot}`` -- plot marker style (default: ``circle``)
* ``--alpha`` -- point opacity (default: 0.8)
* ``--point_size`` -- point size (default: 45)
* ``--legend_location`` -- legend position (default: ``best``)
* ``--format {png,pdf,both}`` -- plot format (default: ``pdf``)
* ``--dpi`` -- PNG resolution (default: 300)
* ``--width`` / ``--height`` -- plot size in inches (default: 8 x 8)
* ``--no_plot`` -- write coordinates without generating a plot
* ``-o``, ``--out_prefix``, ``--output`` -- output prefix

Display all options with::

   beta_tSNE -h


Output
------

For output prefix ``HCV_vs_normal``, the command always writes:

* ``HCV_vs_normal.tSNE.tsv`` -- t-SNE coordinates and group labels

Unless ``--no_plot`` is used, it also writes one or more plots:

* ``HCV_vs_normal.tSNE.pdf``
* ``HCV_vs_normal.tSNE.png``

Use ``--format both`` to generate both plot formats.

The coordinate table contains ``tSNE1`` and ``tSNE2`` for a 2D embedding,
and also ``tSNE3`` when ``--n_components 3`` is used. Plotting always uses the
first two dimensions.


Example Data
------------

* `Beta matrix <https://sourceforge.net/projects/cpgtools/files/test/cirrHCV_vs_normal.data.tsv>`_
* `Group file <https://sourceforge.net/projects/cpgtools/files/test/cirrHCV_vs_normal.grp.csv>`_


Example Figure
--------------

.. image:: ../_static/HCV_vs_normal.tSNE.png
   :height: 450px
   :width: 450px
   :scale: 100%
   :alt: t-SNE plot of DNA methylation samples
