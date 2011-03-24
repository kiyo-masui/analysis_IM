#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


rfi4096_module = Extension('_rfi4096',
                           sources=['rfi4096_wrap.c', 'rfi4096.c'],
                           )

setup (name = 'rfi4096',
       ext_modules = [rfi4096_module],
       py_modules = ["rfi4096"],
       )


