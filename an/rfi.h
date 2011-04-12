double absl(double);
double sq(double);
int get_rms(double *, int, double *, double *);
void xi_square_fit(double *, double *, double *, int);
void get_fit(int, double *, double *, double *);
void flatten(double *, double *, int);
void mask_array(int, double *, double, int *, int);
void rfi_find_dTdf(double *, double *, int *, int, int);
void hi_f_spikes(int, double *, double *, int, int, int);
void rfi_flag(int, int, int, int, double, double *, double *, double *, int *);
void clean(int, double, int, int, int, int, int, double *, double *, 
	       double *, double *, int *);
