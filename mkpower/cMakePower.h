#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926

//#define Linearkbin

typedef struct{
	int dim;

	double *map;
	double *r;	int rn; 
	double *mapinf;
	double ra0;	double dra;
	double dec0;	double ddec;
	int *mapshape;

	double *box;
	int *boxinf0;	int *boxinf1;
	int *boxshape;
	double *box_x;	int nbox_x;
	double *box_y;	int nbox_y;
	double *box_z;	int nbox_z;
}FillConf;


//typedef struct{
//	int dim;
//	double *r;	int rn;
//	double *ra;	int ran;
//	double *dec;int decn;
//	double *box;
//	double *boxinf;
//	int *boxshape;
//	double *map;
//
//	int *mapshape;
//	double *mapinf;
//}FillConf;

typedef struct{
	int dim;
	double *r;	int rn;
	double *ra;	int ran;
	double *dec;int decn;
	double *box;double *box2;
	double *boxinf;
	int *boxshape;
	double *map;double *map2;

	int *mapshape;
	double *mapinf;
}FillConf2;

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

int fillingf2(FillConf2 *conf );

int nfillingf(FillConf2 *conf );

int makepk(FFT *fft, PK *pk);
