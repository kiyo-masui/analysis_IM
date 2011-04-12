%module rfi
%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

%rename(get_fit) get_fit_py;
%rename(clean) clean_py;

int get_fit_py(int DIM1, double* IN_ARRAY1, 
                int DIM1, double* IN_ARRAY1, 
                int DIM1, double* INPLACE_ARRAY1);
int clean_py(double, int, int, int, int, int, int DIM1, double* IN_ARRAY1, 
              int DIM1, double* IN_ARRAY1, int DIM1, double* IN_ARRAY1, 
              int DIM1, double* IN_ARRAY1, int DIM1, int* INPLACE_ARRAY1);

