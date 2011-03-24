
%module rfi4096
%{
#define SWIG_FILE_WITH_INIT
	%}
%include "numpy.i"
%init %{
	import_array();
	%}
void get_fit4096(double INPLACE_ARRAY1[4096], double INPLACE_ARRAY1[4096], double INPLACE_ARRAY1[4096]);
void clean4096(double, int, int, int, int, int, double INPLACE_ARRAY1[4096], double INPLACE_ARRAY1[4096], double INPLACE_ARRAY1[4096], double INPLACE_ARRAY1[4096], int INPLACE_ARRAY1[4096]);

