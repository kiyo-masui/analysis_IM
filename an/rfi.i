%module rfi

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
import_array();
%}

void get_fit(double INPLACE_ARRAY1[2048], double INPLACE_ARRAY1[2048], double INPLACE_ARRAY1[2048]);
void clean(double, int, int, int, int, int, double INPLACE_ARRAY1[2048], double INPLACE_ARRAY1[2048], double INPLACE_ARRAY1[2048], double INPLACE_ARRAY1[2048], int INPLACE_ARRAY1[2048]);
