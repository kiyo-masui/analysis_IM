from distutils.core import setup  
from distutils.extension import Extension  
from Cython.Distutils import build_ext  

import numpy as np

setup(  
   name = 'CubicSpline',  
   ext_modules=[ Extension('cubicspline', ['cubicspline.pyx'], include_dirs=[np.get_include()]) ],  
   cmdclass = {'build_ext': build_ext}  
)
