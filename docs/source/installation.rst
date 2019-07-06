Installation
============================

CpGtools are written in Python. Python3 (v3.5.x) is required to run all programs in
CpGtools. Some programs also need R and R libraries to generate graphs and fit linear and
beta-binomial models.

Prerequisites
--------------
Note: You need to install these tools if they are not available from your computer. 

- `Python 3 <https://www.python.org/downloads/>`_
- `pip3 <https://pip.pypa.io/en/stable/installing/>`_
- `R <https://www.r-project.org/>`_
- R library `aod <https://cran.r-project.org/package-aod>`_ (only required by **dmc_bb.py**)
- R library `beanplot <https://cran.r-project.org/web/packages/beanplot/index.html>`_
  (only needed by **beta_jitter_plot.py**)

Python Dependencies
--------------------
Note: You do NOT need to install these packages manually, as they will be automatically
installed if you use `pip3 <https://pip.pypa.io/en/stable/installing/>`_ to install CpGtools.

- `pandas <https://pandas.pydata.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <https://www.scipy.org/>`_
- `sklearn <https://www.scilearn.com/>`_
- `weblogo <https://pypi.org/project/weblogo/>`_
- `bx-python <https://github.com/bxlab/bx-python>`_

Install CpGtools using pip3 from `PyPI <https://pypi.org/project/cpgtools/>`_ or `github <https://github.com/liguowang/cpgtools>`_
------------------------------------------------------------------------------------------------------------------------------------
::

 $ pip3 install cpgtools
 or 
 $ pip3 install git+https://github.com/liguowang/cpgtools.git
 
Install CpGtools from source code
-----------------------------------
First, download the latest `CpGtools <https://sourceforge.net/projects/cpgtools/files/>`_,
and then execute the following commands

::

 $ tar zxf cpgtools-VERSION.tar.gz
 $ cd cpgtools-VERSION
 $ python3 setup install  #install CpGtools to the default location
 or 
 $ python setup.py install --root-/home/my_pylib/  #install CpGtools to user specified location

After the installation is completed, you probably need to setup up the environment variables
(Below is only an example. Change according to your system configuration)
::

 $ export PYTHONPATH-/home/my_pylib/python3.7/site-packages:$PYTHONPATH

Upgrade CpGtools
-----------------
::

 $ pip3 install cpgtools --upgrade	
