#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
from distutils.core import setup, Extension

import numpy as np

inc = np.get_include()

rfi_module = Extension('_rfi', sources=['rfi_wrap.c', 'rfi.c'], 
                       include_dirs=[inc])

setup (name = 'rfi',
       ext_modules = [rfi_module],
       py_modules = ["rfi"])

