from distutils.core import setup  
from distutils.extension import Extension  
from Cython.Distutils import build_ext  

import numpy as np

setup(  
   name = 'CubicSpline',  
   ext_modules=[ Extension('cubicspline', ['cubicspline.pyx'], include_dirs=[np.get_include(), "/scinet/gpc/Libraries/gsl-1.15-gcc-4.6.1/include"]),
                 Extension('_sphbessel_c', ['_sphbessel_c.pyx'], include_dirs=[np.get_include(), "/scinet/gpc/Libraries/gsl-1.15-gcc-4.6.1/include/"], library_dirs=['/scinet/gpc/Libraries/gsl-1.15-gcc-4.6.1/lib'], libraries=['gsl', 'gslcblas'])],  
   cmdclass = {'build_ext': build_ext}  
)

#old directory was "/opt/gsl-1.15-intel-13.1.3/include"
#old libary_dirs "/opt/gsl-1.15-intel-13.1.3/lib"
