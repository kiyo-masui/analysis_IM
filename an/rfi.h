
extern double absl(double);
extern double sq(double);
extern void testing(double A);
extern void testing2(double A, double B[10]);
extern void get_rms(double A[2048], int, double B[2048], double C[2048]);
extern void xi_square_fit(double A[2048], double B[2048], double C[2048], int);
extern void get_fit(double A[2048], double B[2048], double C[2048]);
extern void flatten(double A[2048], double B[2048]);
extern void mask_array(int, double A[2048], double, int B[2048], int);
extern void rfi_find_dTdf(double A[2048], double B[2048], int C[2048], int, int);
extern void hi_f_spikes(double A[2048], double B[2048], int, int, int);
extern void rfi_flag(int, int, int, double, double A[2048], double B[2048], double C[2048], int D[2048]);
extern void clean(double, int, int, int, int, int, double A[2048], double B[2048], double C[2048], double D[2048], int E[2048]);
