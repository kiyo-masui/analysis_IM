import subprocess

from distutils.core import setup  
from distutils.extension import Extension  
from Cython.Distutils import build_ext  

import numpy as np

cmd = ("gfortran -fopenmp -c -fdefault-real-8 -g cholesky.f90"
       " -o cholesky.o -fPIC")
print cmd
subprocess.check_call(cmd, shell=True)

setup(name='mapmaker',  
     ext_modules=[
                  Extension('chol', ['chol.pyx'], 
                             include_dirs=[np.get_include()], 
                             libraries=[
                                        #"mkl_gnu_thread",  "mkl_core",
                                        #"mkl_intel_ilp64",
                                        "mkl_rt",
                                        'gomp', 'gfortran'],
                             library_dirs=["/scinet/gpc/intel/ics/composer_xe_2015.1.133/mkl/lib/intel64/"],
                             extra_objects=["cholesky.o"],
                             depends=["cholesky.f90", "setup.py",],
                             extra_compile_args=['-fopenmp'],
                             extra_link_args=[],
                             )
                 ],
     cmdclass={'build_ext': build_ext}  
)
