import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for CpGtools  -- Tools to analyze and visualize DNA methylation data (WGBS, RRBS, 850K/450K array).
"""

def main():
    setup(  name = "cpgtools",
            version = "1.0.0",
            python_requires='>=3.5',
            py_modules = [ 'psyco_full' ],
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            package_data = { '': ['*.ps'] },
            scripts = glob.glob( "bin/*.py"),
            ext_modules = [],
            test_suite = 'nose.collector',
            setup_requires = ['nose>=0.10.4'],
            author = "Liguo Wang",
            author_email ="wangliguo78@gmail.com",
            platforms = ['Linux','MacOS'],
            requires = [],
            install_requires = ['numpy','scipy','sklearn','weblogo','bx-python', 'pandas'], 
            description = "Tools to analyze and visualize DNA methylation data",
            url = "http://cpgtools.sourceforge.net/",
            zip_safe = False,
            dependency_links = [],
            classifiers=[
                'Development Status :: 5 - Production/Stable',
                'Environment :: Console',
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
                'Operating System :: MacOS :: MacOS X',
                'Operating System :: POSIX',
                'Programming Language :: Python',
                'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
            
            keywords='DNA methylation, CpG, bisulfite sequencing, RRBS, 850K, 450K',
             )


if __name__ == "__main__":
    main()
