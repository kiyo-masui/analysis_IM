
extern double absl(double);
extern double sq(double);
extern void get_rms(double A[4096], int, double B[4096], double C[4096]);
extern void xi_square_fit(double A[4096], double B[4096], double C[4096], int);
extern void get_fit4096(double A[4096], double B[4096], double C[4096]);
extern void flatten(double A[4096], double B[4096]);
extern void mask_array(int, double A[4096], double, int B[4096], int);
extern void rfi_find_dTdf(double A[4096], double B[4096], int C[4096], int, int);
extern void hi_f_spikes(double A[4096], double B[4096], int, int, int);
extern void rfi_flag(int, int, int, double, double A[4096], double B[4096], double C[4096], int D[4096]);
extern void clean4096(double, int, int, int, int, int, double A[4096], double B[4096], double C[4096], double D[4096], int E[4096]);
