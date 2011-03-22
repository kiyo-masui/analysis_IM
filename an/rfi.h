
extern double absl(double);
extern double sq(double);
extern int get_rms(double *, int, double *, double *);
extern void xi_square_fit(double *, double *, double *, int);
extern void get_fit(int, double *, double *, double *);
extern void flatten(double *, double *, int);
extern void mask_array(int, double *, double, int *, int);
extern void rfi_find_dTdf(double *, double *, int *, int, int);
extern void hi_f_spikes(int, double *, double *, int, int, int);
extern void rfi_flag(int, int, int, int, double, double *, double *, double *, int *);
extern void clean(int, double, int, int, int, int, int, double *, double *, double *, double *, int *);

