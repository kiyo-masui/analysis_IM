import subprocess

from distutils.core import setup  
from distutils.extension import Extension  
from Cython.Distutils import build_ext  

import numpy as np

cmd = ("gfortran -mcmodel=medium -c -fdefault-real-8 -g -O3 cholesky.f90"
       " -o cholesky.o -fPIC")
print cmd
subprocess.check_call(cmd, shell=True)

setup(name='mapmaker',  
     ext_modules=[Extension('_mapmaker', ['_mapmaker.pyx'], 
                            include_dirs=[np.get_include()], 
                            depends=["setup.py",
                                     "_mapmaker.pyx"]),
                  Extension('_cholesky', ['_cholesky.pyx'], 
                             include_dirs=[np.get_include()], 
                             libraries=["pthread", "dl", 
                                        "gfortran"],
                             #libraries=["ml", "pthread", "dl", "CALBLAS",
                             #           "aticalrt", "aticalcl", "gfortran"],
                             library_dirs=["/usr/X11R6/lib"],
                             extra_objects=["cholesky.o"],
                             depends=["cholesky.f90", "setup.py",])
                 ],
     cmdclass={'build_ext': build_ext}  
)
