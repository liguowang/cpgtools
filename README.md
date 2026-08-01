[![PyPI](https://img.shields.io/pypi/v/cpgtools.svg)](https://pypi.org/project/cpgtools/)
[![Python](https://img.shields.io/pypi/pyversions/cpgtools.svg)](https://pypi.org/project/cpgtools/)
[![Documentation](https://readthedocs.org/projects/cpgtools/badge/?version=latest)](https://cpgtools.readthedocs.io/)
[![Downloads](https://static.pepy.tech/badge/cpgtools)](https://pepy.tech/project/cpgtools)

# CpGtools

CpGtools is a Python package for analyzing DNA methylation data generated
from Illumina Infinium arrays and sequencing-based methylation assays. It
provides tools for data preprocessing, quality control, imputation,
epigenetic clock estimation, deconvolution, and downstream analysis.

## Features

- Read and preprocess DNA methylation matrices
- Missing-value imputation
- DNA methylation age estimation (epigenetic clocks)
- Cell-type deconvolution
- Differential methylation analysis
- Data visualization

## Installation

Install the latest stable release from PyPI:

```bash
pip install cpgtools
```

Or install the latest development version from GitHub:

```bash
pip install git+https://github.com/liguowang/cpgtools.git
```

### (Optional) Create a virtual environment

Using a virtual environment is recommended but not required.

```bash
python3 -m venv cpgtools-env
source cpgtools-env/bin/activate
```

## Upgrade

```bash
pip install --upgrade cpgtools
```

## Uninstall

```bash
pip uninstall -y cpgtools
```

## Documentation

Documentation is available at:

https://cpgtools.readthedocs.io/

## Source Code

https://github.com/liguowang/cpgtools

## License

CpGtools is distributed under the MIT License.
