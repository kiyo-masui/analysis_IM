#include <math.h>
#include <fftw.h>
#include "cFFTW.h"

int fft3d(fftDATA *ff){
	fftw_complex * data =
		(fftw_complex *)malloc(
			ff->sizex * ff->sizey * ff->sizez * sizeof(fftw_complex));
	for(int i=0; i<(ff->sizex * ff->sizey * ff->sizez); i++){
		data[i].re = ff->data[i];
		data[i].im = 0.;
	}
	
	fftwnd_plan p = fftw3d_create_plan(ff->sizex, ff->sizey, ff->sizez, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE);

	fftwnd_one(p, data, NULL);

	for(int i=0; i<(ff->sizex * ff->sizey * ff->sizez); i++){
		double re = data[i].re;
		double im = data[i].im;
		ff->data[i] = sqrt(re*re + im*im);
//		if(!(ff->data[i][j][k]>1.||ff->data[i][j][k]<=1.))
//			printf("%lg %lg\n",re, im);
	}

	fftwnd_destroy_plan(p);

	return 0;

}

//int fft3dc(fftDATA *ff, fftDATA *ff2){
//	fftwnd_plan p = fftw3d_create_plan(ff->size, ff->size, ff->size, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE);
//
//	fftw_complex * data =
//		(fftw_complex *)malloc(pow(ff->size, ff->dim)*sizeof(fftw_complex));
//	fftw_complex * data2 =
//		(fftw_complex *)malloc(pow(ff->size, ff->dim)*sizeof(fftw_complex));
//
//	for(int i=0; i<pow(ff->size, ff->dim); i++){
//		data[i].re = ff->data[i];
//		data2[i].re = ff2->data[i];
//		data[i].im = 0.;
//		data2[i].im = 0.;
//	}
//	
//
//	fftwnd_one(p, data, NULL);
//	fftwnd_one(p, data2, NULL);
//
//	for(int i=0; i<pow(ff->size, ff->dim); i++){
//		double re = data[i].re;
//		double re2 = data2[i].re;
//		double im = data[i].im;
//		double im2 = data2[i].im;
//		ff->data[i] = 
//			sqrt(pow((re*re2+im*im2),2) + pow((re*im2-re2*im),2));
////		if(!(ff->data[i][j][k]>1.||ff->data[i][j][k]<=1.))
////			printf("%lg %lg\n",re2, im2);
//	}
//
//	fftwnd_destroy_plan(p);
//
//	return 0;
//
//}
