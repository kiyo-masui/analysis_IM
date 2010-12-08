#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


rfi_module = Extension('_rfi',
                           sources=['rfi_wrap.c', 'rfi.c'],
                           )

setup (name = 'rfi',
       ext_modules = [rfi_module],
       py_modules = ["rfi"],
       )


