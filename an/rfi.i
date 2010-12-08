%module rfi
%include "carrays.i"
%array_functions(double, doublearray);
%array_functions(int, intarray);

%{
#define SWIG_FILE_WITH_INIT
%}

double absl(double);
double sq(double);
void testing(double A);
void testing2(double A, double B[10]);
void get_rms(double A[2048], int, double B[2048], double C[2048]);
void xi_square_fit(double A[2048], double B[2048], double C[2048], int);
void get_fit(double A[2048], double B[2048], double C[2048]);
void flatten(double A[2048], double B[2048]);
void mask_array(int, double A[2048], double, int B[2048], int);
void rfi_find_dTdf(double A[2048], double B[2048], int C[2048], int, int);
void hi_f_spikes(double A[2048], double B[2048], int, int, int);
void rfi_flag(int, int, int, double, double A[2048], double B[2048], double C[2048], int D[2048]);
void clean(double, int, int, int, int, int, double A[2048], double B[2048], double C[2048], double D[2048], int E[2048]);
