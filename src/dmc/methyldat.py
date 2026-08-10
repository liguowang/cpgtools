#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 10:55:14 2022
"""
from dmc import ireader
import os

class MakeMethylObj():
    """
    Pack DNA methylation signature file into object.

    >>> from dmc.methyldat import MakeMethylObj
    >>> a = MakeMethylObj(signature_file = 'coefBlup.tsv', signature_name = \
                          'BLUP', signature_info="BLUP clock", unit='year')
    >>> a.name
    'BLUP'
    >>> a.Intercept
    91.15396
    >>> a.ncpg
    319607
    >>> a.cpgs[1:5]
    ['cg14361672', 'cg01763666', 'cg02115394', 'cg13417420']
    >>> a.coef['cg14361672']
    0.000557006779951299
    """

    def __init__(self, signature_file, signature_name, tissues=[],
                 unit='', signature_info='', reference='', pub_link='',
                 method='', organism='', training_platform=[],
                 training_age_range=[], training_age_range_unit=''):
        self.name = signature_name
        self.info = signature_info
        self.tissues = tissues
        self.unit = unit
        self.coef = {}
        self.cpgs = []
        self.ncpg = 0
        self.Intercept = 0.0
        self.ref = reference
        self.pubmed = pub_link
        self.method = method
        self.organism = organism
        self.t_platform = training_platform
        self.age_range = training_age_range
        self.age_unit = training_age_range_unit
        if os.path.exists(signature_file) and os.path.getsize(signature_file) > 0:
            for l in ireader.reader(signature_file):
                if l.startswith('#'):
                    continue
                f = l.split()
                if l.startswith('Intercept'):
                    try:
                        self.Intercept = float(f[1])
                    except:
                        self.Intercept = 0.0
                else:
                    self.cpgs.append(f[0])
                    self.ncpg  += 1
                    try:
                        self.coef[f[0]] = float(f[1])
                        # self.ncpg  += 1
                    except:
                        continue
        else:
            self.Intercept = 'N/A'
            self.cpgs = []
            self.ncpg = 'N/A'
            self.coef = []
