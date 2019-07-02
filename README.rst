Install CpGtools using pip3 
----------------------------
::

 $ pip3 install cpgtools
 or 
 $ pip3 install git+https://github.com/liguowang/cpgtools.git
 
Instal CpGtools from source code
--------------------------------
First download the latest `CpGtools <https://sourceforge.net/projects/cpgtools/files/>`_,
and then execute the following commands

::

 $ tar zxf cpgtools-VERSION.tar.gz
 $ cd cpgtools-VERSION
 $ python3 setup install	#install CpGtools to the default location
 or 
 $ python setup.py install --root=/home/my_pylib/	#install CpGtools to user specified location

After the installation is completed, you probably need to setup up the environment variables
(Below is only an example. Change according to your system configuration)
::

 $ export PYTHONPATH=/home/my_pylib/python3.7/site-packages:$PYTHONPATH

Upgrade CpGtools
-----------------
::

 $ pip3 install cpgtools --upgrade	


Online manual
--------------
https://cpgtools.readthedocs.io/en/latest/
