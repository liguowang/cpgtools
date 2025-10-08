Installation
============

**CpGtools** is written in Python and requires **Python 3.5** or later.  
Some tools within CpGtools also require **R** and specific **R libraries** to generate plots and fit linear or beta-binomial models.

Prerequisites
--------------

If these dependencies are not already installed on your system, please install them first:

- `Python 3 <https://www.python.org/downloads/>`_
- `pip3 <https://pip.pypa.io/en/stable/installing/>`_
- `R <https://www.r-project.org/>`_
- R package `aod <https://cran.r-project.org/package=aod>`_ (required by **dmc_bb.py**)
- R package `beanplot <https://cran.r-project.org/web/packages/beanplot/index.html>`_ (required by **beta_jitter_plot.py**)

Dependencies
-------------------

You do **not** need to install these manually — `pip3` will automatically install them when you install CpGtools.

- `pandas <https://pandas.pydata.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <https://www.scipy.org/>`_
- `scikit-learn <https://scikit-learn.org/>`_
- `weblogo <https://pypi.org/project/weblogo/>`_
- `umap-learn <https://pypi.org/project/umap-learn/>`_
- `fancyimpute <https://github.com/iskandr/fancyimpute>`_
- `bx-python <https://github.com/bxlab/bx-python>`_
- `pycombat <https://github.com/epigenelabs/pyComBat>`_
- `matplotlib <https://matplotlib.org/>`_

Install CpGtools using pip3
---------------------------

You can install CpGtools directly from **PyPI** or **GitHub**:

.. code-block:: bash

   $ pip3 install cpgtools
   # or
   $ pip3 install git+https://github.com/liguowang/cpgtools.git

Install CpGtools from Source
----------------------------

To install CpGtools from source:

1. Download the latest release from  
   `CpGtools on SourceForge <https://sourceforge.net/projects/cpgtools/files/>`_.
2. Extract and install:

.. code-block:: bash

   $ tar zxf cpgtools-VERSION.tar.gz
   $ cd cpgtools-VERSION
   $ python3 setup.py install          # install to the default location
   # or
   $ python3 setup.py install --root=/home/my_pylib/  # install to a custom location

After installation, you may need to set up your environment variables.  
Below is an example (adjust the path to match your setup):

.. code-block:: bash

   $ export PYTHONPATH=/home/my_pylib/python3.7/site-packages:$PYTHONPATH

Upgrade CpGtools
----------------

To upgrade to the latest version:

.. code-block:: bash

   $ pip3 install --upgrade cpgtools

