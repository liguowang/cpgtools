epical
======

Overview
--------

``epical`` provides a unified command-line interface to DNA methylation
epigenetic clocks and related methylation-based predictors distributed with
CpGtools.

The command uses a **subcommand-per-clock** design. Each supported clock has
its own help page and can be run as::

   epical CLOCK Input_file [options]

For example::

   epical Horvath13 beta.tsv -o sample

Display the top-level help with::

   epical -h

Display help for a specific clock with::

   epical Horvath13 -h


Supported Clocks
----------------

The current CLI exposes the following clocks and predictors.

Human and general clocks
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Command
     - Notes
   * - ``Horvath13``
     - Standard epigenetic clock command.
   * - ``Horvath13_shrunk``
     - Shrunk version of the Horvath 2013 clock.
   * - ``Horvath18``
     - Horvath 2018 clock.
   * - ``Levine``
     - Levine methylation clock.
   * - ``Hannum``
     - Hannum methylation clock.
   * - ``Zhang_EN``
     - Zhang elastic-net clock.
   * - ``Zhang_BLUP``
     - Zhang BLUP clock.
   * - ``AltumAge``
     - AltumAge predictor.
   * - ``Lu_DNAmTL``
     - DNAm telomere-length predictor.
   * - ``Weidner``
     - Weidner clock.
   * - ``Lin``
     - Lin clock.
   * - ``ENCen100``
     - EN-Cen 100-CpG model.
   * - ``ENCen40``
     - EN-Cen 40-CpG model.
   * - ``Ped_Wu``
     - Pediatric Wu clock.
   * - ``PedBE``
     - Pediatric buccal epigenetic clock.
   * - ``Cortical``
     - Cortical methylation clock.
   * - ``MEAT``
     - Muscle epigenetic age test.
   * - ``mammClock1``
     - Mammalian clock using the general mammalian implementation.

Gestational-age clocks
~~~~~~~~~~~~~~~~~~~~~~

The following gestational-age predictors are available:

* ``GA_Bohlin``
* ``GA_Haftorn``
* ``GA_Knight``
* ``GA_Mayne``
* ``GA_Lee_CPC``
* ``GA_Lee_RPC``
* ``GA_Lee_rRPC``


DunedinPACE
~~~~~~~~~~~

``DunedinPACE`` is available as a dedicated subcommand. Its CLI is similar to
the standard clocks, except that it intentionally does **not** provide the
``--log`` option.


Mouse clocks
~~~~~~~~~~~~

The following mouse clocks are available:

* ``WLMT``
* ``YOMT``
* ``mmLiver``
* ``mmBlood``

These commands add a genome-build option::

   -g {mm10,mm39}
   --genome {mm10,mm39}

The default is ``mm10``.


Mammalian species-specific clocks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following clocks accept an explicit mammalian species:

* ``mammClock2``
* ``mammClock3``

Use::

   -s {human,mouse}
   --species {human,mouse}

The default is ``human``.


EPM
~~~

``EPM`` uses a separate model-building interface and requires both a
methylation matrix and a metadata file as positional arguments::

   epical EPM Input_file meta_file [options]

See :ref:`epm-options` below for its specific parameters.


Input
-----

For most clock subcommands, the first positional argument is the methylation
input file::

   epical Horvath13 Input_file

The exact accepted matrix format and delimiter behavior are defined by the
underlying CpGtools methylation-clock implementation. Use the command-specific
help for details::

   epical Horvath13 -h


Common Options
--------------

Most clock commands share the following options:

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Option
     - Description
   * - ``-o``, ``--output``
     - Output prefix.
   * - ``-p``, ``--percent``
     - Missing-CpG threshold used by the clock implementation.
       Default: ``0.2``.
   * - ``-d``, ``--delimiter``
     - Input delimiter. When omitted, the underlying reader determines the
       delimiter.
   * - ``-f {pdf,png}``, ``--format {pdf,png}``
     - Figure output format. Default: ``pdf``.
   * - ``-m``, ``--metadata``
     - Optional metadata file.
   * - ``-l``, ``--log``
     - Optional log-file path. This option is not available for
       ``DunedinPACE``.
   * - ``--impute``
     - Imputation method code. Accepted values are integers from ``-1`` through
       ``11``. Default: ``11``.
   * - ``-r``, ``--ref``
     - Optional external reference file used by supported imputation methods.
   * - ``--overwrite``
     - Allow existing output files to be overwritten.
   * - ``--debug``
     - Enable debug logging.


Example
-------

Run the Horvath 2013 clock::

   epical Horvath13 beta.tsv \
       -o horvath13

Specify a metadata file and PNG output::

   epical Horvath13 beta.tsv \
       -m metadata.tsv \
       -f png \
       -o horvath13

Use an external reference file for imputation::

   epical Horvath13 beta.tsv \
       --impute 11 \
       -r reference.tsv \
       -o horvath13


Mouse-clock Example
-------------------

Run a mouse clock using the mm39 genome build::

   epical WLMT mouse_beta.tsv \
       --genome mm39 \
       -o mouse_age


Mammalian-clock Example
-----------------------

Run a mammalian clock for mouse samples::

   epical mammClock2 beta.tsv \
       --species mouse \
       -o mammalian_age


.. _epm-options:

EPM Options
-----------

``EPM`` has its own interface::

   epical EPM Input_file meta_file [options]

Its options are:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``-o``, ``--output``
     - Output prefix.
   * - ``-p``, ``--pcc``
     - Absolute Pearson correlation coefficient threshold between
       chronological age and Beta-values. Default: ``0.85``.
   * - ``-n``, ``--niter``
     - Number of expectation-maximization iterations. Default: ``100``.
   * - ``-k``, ``--kfold``
     - Number of cross-validation folds. Default: ``10``.
   * - ``-e``, ``--etol``
     - Model-fitting error tolerance. Default: ``1e-5``.
   * - ``-d``, ``--delimiter``
     - Input delimiter.
   * - ``-f {pdf,png}``, ``--format {pdf,png}``
     - Figure output format. Default: ``pdf``.
   * - ``-l``, ``--log``
     - Optional log-file path.
   * - ``--impute``
     - Imputation method code from ``-1`` through ``11``. Default: ``11``.
   * - ``-r``, ``--ref``
     - Optional external reference file.
   * - ``--debug``
     - Enable debug logging.

Example::

   epical EPM beta.tsv metadata.tsv \
       --pcc 0.85 \
       --niter 100 \
       --kfold 10 \
       -o epm_model


Version
-------

Display the installed version with::

   epical --version


Command-specific Help
---------------------

Each clock has its own help text derived from the model metadata bundled with
CpGtools. Because clock-specific assumptions and requirements may differ,
check the relevant subcommand before running an analysis.

Examples::

   epical Horvath13 -h
   epical DunedinPACE -h
   epical WLMT -h
   epical mammClock2 -h
   epical EPM -h
