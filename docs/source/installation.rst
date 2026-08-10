Installation
============

**CpGtools** is a Python package for analyzing DNA methylation data. It runs on
Linux, macOS, and Windows, and requires **Python 3.9 or later**.

Most Python dependencies are installed automatically when CpGtools is
installed with ``pip``. Some commands additionally require **R** and specific
R packages.

Prerequisites
-------------

Before installing CpGtools, ensure that the following software is available on
your system:

* `Python 3 <https://www.python.org/downloads/>`_
* `pip <https://pip.pypa.io/en/stable/installation/>`_

The following software is optional and is required only by specific
subcommands:

* `TensorFlow <https://www.tensorflow.org/>`_
  (required by ``beta_impute.py morel --model DNN``)
* `R <https://www.r-project.org/>`_
* R package `aod <https://cran.r-project.org/package=aod>`_
  (required by ``dmc_bb.py``)
* R package
  `beanplot <https://cran.r-project.org/package=beanplot>`_
  (required by ``beta_jitter_plot.py``)



Install from PyPI
-----------------

The recommended way to install CpGtools is from the Python Package Index
(PyPI):

.. code-block:: bash

   pip install cpgtools

Install the Development Version
-------------------------------

To install the latest development version directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/liguowang/cpgtools.git

Install from Source
-------------------

Clone the repository and install CpGtools locally:

.. code-block:: bash

   git clone https://github.com/liguowang/cpgtools.git
   cd cpgtools
   pip install .

Alternatively, install in editable (development) mode:

.. code-block:: bash

   pip install -e .

Dependencies
------------

The following Python packages are installed automatically when CpGtools is
installed via ``pip``:

* ``numpy``
* ``scipy``
* ``pandas``
* ``matplotlib``
* ``scikit-learn``
* ``fancyimpute``
* ``umap-learn``
* ``weblogo``
* ``bx-python``
* ``pycombat``


Optional Python Dependencies
----------------------------

Most Python dependencies are installed automatically when CpGtools is
installed via ``pip``. The following package is optional:

* ``tensorflow`` — required only when running the MOREL algorithm with
  ``--model DNN`` in ``beta_impute.py``. It is **not** required for the
  default Random Forest implementation (``--model RF``).

To install TensorFlow:

.. code-block:: bash

   pip install tensorflow


Upgrade CpGtools
----------------

To upgrade an existing installation to the latest release:

.. code-block:: bash

   pip install --upgrade cpgtools

Verify the Installation
-----------------------

To verify that CpGtools was installed successfully:

.. code-block:: bash

   cpgtools --version

or

.. code-block:: bash

   python -m cpgtools --version
