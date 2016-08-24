import subprocess

from distutils.core import setup  
from distutils.extension import Extension  
from Cython.Distutils import build_ext  

import numpy as np

#cmd = ("gfortran -mcmodel=medium -c -fdefault-real-8 -fdefault-integer-8  -g -O3 cholesky.f90"
#       " -o cholesky.o -fPIC -m64")

#cmd = ("gfortran -mcmodel=medium -c -fdefault-real-8 -fdefault-integer-8  -ffree-form -finteger-4-integer-8 -w -fbounds-check  -g3 -O0 -fbacktrace -Wuninitialized cholesky.f90"
#       " -o cholesky.o -fPIC -m64")

#cmd = ("ifort -i8 -r8 -mkl cholesky.f90" " -o cholesky.o -fPIC -m64" )

#cmd = ("gfortran -fopenmp -fdefault-real-8 -fdefault-integer-8  -c cholesky.f90 -L/scinet/gpc/intel/ics/composer_xe_2015.1.133/mkl/lib/intel64/ -lmkl_gnu_thread  -lmkl_core  -lmkl_intel_ilp64" " -o cholesky.o -fPIC -m64")

cmd = ("gfortran -fopenmp -c -fdefault-real-8 -g cholesky.f90"
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
                             #libraries=["mkl_intel_ilp64", "mkl_core",
                             #"mkl_gnu_thread","gfortran","dl","m",
                             #"pthread","gomp"],
                             libraries =["mkl_rt","gfortran","dl","m",
                             "pthread","gomp"], 
                             #libraries=["pthread", "dl", 
                             #           "gfortran"],
                             #libraries=["ml", "pthread", "dl", "CALBLAS",
                             #           "aticalrt", "aticalcl", "gfortran"],
                             extra_compile_args = ["-fopenmp", "-m64"],
                             extra_link_args = ["-Wl,--no-as-needed"],
                             #library_dirs=["/usr/X11R6/lib"],
                             library_dirs=["/scinet/gpc/intel/ics/composer_xe_2015.1.133/mkl/lib/intel64/"],
                             extra_objects=["cholesky.o"],
                             depends=["cholesky.f90", "setup.py",])
                 ],
     cmdclass={'build_ext': build_ext}  
)
