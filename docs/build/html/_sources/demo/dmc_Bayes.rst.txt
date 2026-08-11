dmc_Bayes
=========

Overview
--------

``dmc_Bayes`` estimates methylation differences between two groups using a
Bayesian Markov chain Monte Carlo (MCMC) approach inspired by John Kruschke's
BEST framework (Bayesian Estimation Supersedes the t test).

For each CpG/probe, the command estimates:

* posterior mean for group 1;
* posterior mean for group 2;
* posterior difference of means;
* 95% highest density interval (HDI) for the difference; and
* posterior probability that the difference is on the same side of zero as
  the posterior median.

Because an independent MCMC chain is fitted for every CpG, this method is much
slower than a standard t test. Multiple worker processes are recommended for
large datasets.


Input Files
-----------

Methylation matrix
~~~~~~~~~~~~~~~~~~

The input matrix must have **CpGs/probes in rows** and **samples in columns**.
The first row contains sample IDs and the first column contains CpG/probe IDs.

Example::

   CpG_ID      Sample_01   Sample_02   Sample_03   Sample_04
   cg00001099  0.73        0.81        0.79        0.84
   cg00000363  0.62        0.59        0.45        0.49

Non-numeric and non-finite values are treated as missing and ignored.
Plain-text, ``.gz``, and ``.bz2`` inputs are supported.


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

All samples listed in the group file must be present in the methylation
matrix. Samples present in the matrix but absent from the group file are not
used.

The two group IDs are sorted internally; ``mu1`` corresponds to the first
group in sorted order and ``mu2`` to the second.


Usage
-----

Basic usage::

   dmc_Bayes \
       -i test_05_TwoGroup.tsv.gz \
       -g test_05_TwoGroup.grp.csv \
       -p 10 \
       -o dmc_output

Useful options include:

* ``-n``, ``--niter`` -- number of MCMC iterations (default: 5000)
* ``-b``, ``--burnin`` -- number of initial MCMC iterations discarded as
  burn-in (default: 500)
* ``-p``, ``--processor`` -- number of worker processes (default: 1)
* ``-s``, ``--seed`` -- random-number seed (default: 99)
* ``-o``, ``--output`` -- output prefix

``--niter`` must be greater than 1, ``--burnin`` must be non-negative and
smaller than ``--niter``, and ``--processor`` must be positive.

Display all options with::

   dmc_Bayes -h


Output
------

For output prefix ``dmc_output``, the command writes:

* ``dmc_output.bayes.tsv``

The output contains six columns:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Description
   * - ``ID``
     - CpG/probe ID.
   * - ``mu1``
     - Posterior mean methylation level for group 1.
   * - ``mu2``
     - Posterior mean methylation level for group 2.
   * - ``mu_diff``
     - Posterior mean difference, ``mu1 - mu2``.
   * - ``mu_diff (95% HDI)``
     - 95% highest density interval of the posterior difference.
   * - ``Probability``
     - Posterior probability that the difference has the same sign as its
       posterior median. For a positive median this is
       :math:`P(\mu_1-\mu_2>0)`; for a negative median this is
       :math:`P(\mu_1-\mu_2<0)`.

Example::

   ID           mu1      mu2      mu_diff    mu_diff (95% HDI)      Probability
   cg00001099   0.775209 0.795404 -0.020196  (-0.065148,0.023974)   0.811024
   cg00000363   0.610565 0.469523  0.141042  (0.030769,0.232965)    0.994665


Interpretation
--------------

A 95% HDI that does not include zero indicates that the posterior difference
is credibly separated from zero under this model.

``Probability`` should not be interpreted as a frequentist p-value. It
summarizes the posterior directional probability of the estimated difference.

A CpG with fewer than two valid observations in either group is reported with
``nan`` values because the model cannot be fitted for that CpG.


Reproducibility
---------------

The random seed is deterministic. Each CpG receives a separate seed derived
from ``--seed`` and its row order, so results remain reproducible when using
multiple worker processes.


Example Data
------------

* `Methylation matrix
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
* `Group file
  <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_
