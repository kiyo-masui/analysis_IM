typedef struct{
	double *data;
	int dim;
	int sizex;
	int sizey;
	int sizez;
}fftDATA;

int fft3d(fftDATA *ff);

int fft3dc(fftDATA *ff, fftDATA *ff2);
