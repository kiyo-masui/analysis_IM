#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926

typedef struct{
	int dim;
	double *r;	int rn;
	double *ra;	int ran;
	double *dec;int decn;
	double *box; double *box2; 
	double *boxinf;
	int *boxshape;
	double *map; double *map2;

	int *mapshape;
	double *mapinf;
}FillConf;

typedef struct{
	int dim;
	int sizex;
	int sizey;
	int sizez;
	double *data;
}FFT;

typedef struct{
	int N;
	int Np;
	int Nv;
	double *val;
	double *k;
	double *val2;
	double *k2;
}PK;

int fillingf(FillConf *conf );

int nfillingf(FillConf *conf );

int makepk(FFT *fft, PK *pk);
