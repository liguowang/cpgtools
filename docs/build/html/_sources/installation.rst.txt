Installation
============

**CpGtools** supports Python 3.9 or later and can be installed with ``pip``. Most Python
dependencies are installed automatically.

Some commands require additional software, such as R or TensorFlow.



Install in a Virtual Environment
--------------------------------

Using a virtual environment is recommended because it keeps CpGtools and its
Python dependencies isolated from the system Python installation and from
other projects.

Create a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, create a new environment using Python's built-in ``venv`` module:

.. code-block:: bash

   python3 -m venv cpgtools-env

Activate the environment:

On Linux or macOS:

.. code-block:: bash

   source cpgtools-env/bin/activate

On Windows Command Prompt:

.. code-block:: bat

   cpgtools-env\Scripts\activate.bat

On Windows PowerShell:

.. code-block:: powershell

   cpgtools-env\Scripts\Activate.ps1

After activation, the environment name usually appears at the beginning of
the command prompt, for example::

   (cpgtools-env) $

Install CpGtools
~~~~~~~~~~~~~~~~

Upgrade ``pip`` and install CpGtools from PyPI:

.. code-block:: bash

   python -m pip install --upgrade pip
   python -m pip install cpgtools

Verify the installation:

.. code-block:: bash

   epical --version
   beta_impute --version

When finished, leave the virtual environment with:

.. code-block:: bash

   deactivate

To use CpGtools again later, reactivate the same environment rather than
reinstalling the package.

Install from PyPI
-----------------

The recommended installation method is:

.. code-block:: bash

   python -m pip install cpgtools

To upgrade an existing installation:

.. code-block:: bash

   python -m pip install --upgrade cpgtools


Install from GitHub
-------------------

To install the latest development version directly from GitHub:

.. code-block:: bash

   python -m pip install git+https://github.com/liguowang/cpgtools.git


Install from Source
-------------------

Clone the repository and install CpGtools locally:

.. code-block:: bash

   git clone https://github.com/liguowang/cpgtools.git
   cd cpgtools
   python -m pip install .

For development, use an editable installation:

.. code-block:: bash

   python -m pip install -e .


Python Dependencies
-------------------

The following packages are installed automatically with CpGtools:

* ``numpy``
* ``scipy``
* ``pandas``
* ``scikit-learn``
* ``matplotlib``
* ``umap-learn``
* ``bx-python``
* ``weblogo``
* ``pycombat``

Additional dependencies required by these packages are resolved automatically
by ``pip``.


Optional Dependencies
---------------------

TensorFlow
~~~~~~~~~~

TensorFlow is required only when using the MOREL imputation method with the
dense neural-network model:

.. code-block:: bash

   beta_impute morel --model DNN ...

The default Random Forest model (``--model RF``) does not require TensorFlow.

Install TensorFlow separately if needed:

.. code-block:: bash

   python -m pip install tensorflow


R and R Packages
~~~~~~~~~~~~~~~~

Some CpGtools commands call R and therefore require an R installation.

`R <https://www.r-project.org/>`_
   Required by commands that execute generated R scripts.

`aod <https://cran.r-project.org/package=aod>`_
   Required by ``dmc_bb``.

`beanplot <https://cran.r-project.org/package=beanplot>`_
   Required by ``beta_jitter_plot``.

These R dependencies are not installed automatically by ``pip``.


Verify the Installation
-----------------------

After installation, verify several command-line programs:

.. code-block:: bash

   epical --version
   beta_impute --version
   epical -h
   beta_impute -h
   beta_deconvolution -h

You can also verify that the Python package is importable:

.. code-block:: bash

   python -c "import cpgmodule; print(cpgmodule.__file__)"


Troubleshooting
---------------

If a command is not found after installation, confirm that CpGtools is
installed in the active Python environment:

.. code-block:: bash

   python -m pip show cpgtools

When using a virtual or Conda environment, make sure that environment is
activated before installing or running CpGtools.
