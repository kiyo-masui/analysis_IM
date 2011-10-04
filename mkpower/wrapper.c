#include <Python.h>
#include <arrayobject.h>
#include "cMakePower.h"
//#include "cFFTW.h"

static PyObject * Filling(PyObject *self, PyObject *args){

	FillConf *conf = (FillConf *)malloc(sizeof(FillConf));
	PyArrayObject *map;
	PyArrayObject *r;
	PyArrayObject *mapinf;
	PyArrayObject *box;
	PyArrayObject *boxinf0;
	PyArrayObject *boxinf1;
	PyArrayObject *box_x;
	PyArrayObject *box_y;
	PyArrayObject *box_z;
	
	if(!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!", 
		&PyArray_Type, &map,	
		&PyArray_Type, &r, 
		&PyArray_Type, &mapinf,
		&PyArray_Type, &box, 
		&PyArray_Type, &boxinf0,
		&PyArray_Type, &boxinf1,
		&PyArray_Type, &box_x,
		&PyArray_Type, &box_y,
		&PyArray_Type, &box_z
		))
		return NULL;
	conf->dim = map->nd;
	conf->map = (double *)map->data; // get the map data 
	conf->r = (double *)r->data;		conf->rn = r->dimensions[0];
	conf->mapinf = (double *)mapinf->data;
	conf->ra0 = conf->mapinf[0];
	conf->dra = conf->mapinf[1];
	conf->dec0 = conf->mapinf[2];
	conf->ddec = conf->mapinf[3];
	conf->mapshape = (int *)malloc(conf->dim*sizeof(int));

	conf->box = (double *)box->data; 
	conf->boxshape = (int *)malloc(conf->dim*sizeof(int));
	for (int i=0; i<conf->dim; i++){
		conf->mapshape[i] = map->dimensions[i];
		conf->boxshape[i] = box->dimensions[i];
	}
	conf->boxinf0 = (int *)boxinf0->data;
	conf->boxinf1 = (int *)boxinf1->data;
	conf->box_x = (double *)box_x->data; conf->nbox_x = box_x->dimensions[0];
	conf->box_y = (double *)box_y->data; conf->nbox_y = box_y->dimensions[0];
	conf->box_z = (double *)box_z->data; conf->nbox_z = box_z->dimensions[0];

	fillingf(conf);

//	double MAX = -1.e10; double MIN = 1.e10;
//	double MAX2 = -1.e10; double MIN2 = 1.e10;
//	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
//		if(conf->box[i]>MAX) MAX = conf->box[i];
//		if(conf->box[i]<MIN) MIN = conf->box[i];
//		if(conf->box2[i]>MAX2) MAX2 = conf->box2[i];
//		if(conf->box2[i]<MIN2) MIN2 = conf->box2[i];
//	}
//	printf("MAX = %lg , MIN = %lg\n", MAX, MIN);
//	printf("MAX2 = %lg , MIN2 = %lg\n", MAX2, MIN2);
	return Py_BuildValue("i", 0);
}

static PyObject * Filling2(PyObject *self, PyObject *args){

	FillConf2 *conf = (FillConf2 *)malloc(sizeof(FillConf2));
	PyArrayObject *map;
	PyArrayObject *map2;
	PyArrayObject *box;
	PyArrayObject *box2;
	PyArrayObject *r;
	PyArrayObject *ra;
	PyArrayObject *dec;
	PyArrayObject *boxinf;
	PyArrayObject *mapinf;
	
	if(!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!", 
		&PyArray_Type, &map,	&PyArray_Type, &map2,
		&PyArray_Type, &box, &PyArray_Type, &box2,
		&PyArray_Type, &r, &PyArray_Type, &ra,	&PyArray_Type, &dec,
		&PyArray_Type, &boxinf, &PyArray_Type, &mapinf))
		return NULL;
	conf->dim = map->nd;
	conf->map = (double *)map->data; conf->map2 = (double *)map2->data;
	conf->mapshape = (int *)malloc(conf->dim*sizeof(int));
	conf->boxshape = (int *)malloc(conf->dim*sizeof(int));
	for (int i=0; i<conf->dim; i++){
		conf->mapshape[i] = map->dimensions[i];
		conf->boxshape[i] = box->dimensions[i];
	}
	conf->boxinf = (double *)boxinf->data;
	conf->mapinf = (double *)mapinf->data;
	conf->r = (double *)r->data;		conf->rn = r->dimensions[0];
	conf->ra = (double *)ra->data;	conf->ran = ra->dimensions[0];
	conf->dec = (double *)dec->data;	conf->decn = dec->dimensions[0];
	conf->box = (double *)box->data; conf->box2 = (double *)box2->data;
	fillingf2(conf);

//	double MAX = -1.e10; double MIN = 1.e10;
//	double MAX2 = -1.e10; double MIN2 = 1.e10;
//	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
//		if(conf->box[i]>MAX) MAX = conf->box[i];
//		if(conf->box[i]<MIN) MIN = conf->box[i];
//		if(conf->box2[i]>MAX2) MAX2 = conf->box2[i];
//		if(conf->box2[i]<MIN2) MIN2 = conf->box2[i];
//	}
//	printf("MAX = %lg , MIN = %lg\n", MAX, MIN);
//	printf("MAX2 = %lg , MIN2 = %lg\n", MAX2, MIN2);
	return Py_BuildValue("i", 0);
}

static PyObject * nFilling(PyObject *self, PyObject *args){

	FillConf2 *conf = (FillConf2 *)malloc(sizeof(FillConf2));
	PyArrayObject *map;
	PyArrayObject *map2;
	PyArrayObject *box;
	PyArrayObject *box2;
	PyArrayObject *r;
	PyArrayObject *ra;
	PyArrayObject *dec;
	PyArrayObject *boxinf;
	PyArrayObject *mapinf;
	
	if(!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!", 
		&PyArray_Type, &map,	&PyArray_Type, &map2,
		&PyArray_Type, &box, &PyArray_Type, &box2,
		&PyArray_Type, &r, &PyArray_Type, &ra,	&PyArray_Type, &dec,
		&PyArray_Type, &boxinf, &PyArray_Type, &mapinf))
		return NULL;
	conf->dim = map->nd;
	conf->map = (double *)map->data; conf->map2 = (double *)map2->data;
	conf->box = (double *)box->data; conf->box2 = (double *)box2->data;
	conf->mapshape = (int *)malloc(conf->dim*sizeof(int));
	conf->boxshape = (int *)malloc(conf->dim*sizeof(int));
	for (int i=0; i<conf->dim; i++){
		conf->mapshape[i] = map->dimensions[i];
		conf->boxshape[i] = box->dimensions[i];
	}
	conf->boxinf = (double *)boxinf->data;
	conf->mapinf = (double *)mapinf->data;
	conf->r = (double *)r->data;		conf->rn = r->dimensions[0];
	conf->ra = (double *)ra->data;	conf->ran = ra->dimensions[0];
	conf->dec = (double *)dec->data;	conf->decn = dec->dimensions[0];
	nfillingf(conf);
	return Py_BuildValue("i", 0);
}
static PyObject * Make(PyObject *self, PyObject *args){
	
	FFT *fft = (FFT *)malloc(sizeof(FFT));
	PK *pk = (PK *)malloc(sizeof(PK));
	PyArrayObject *Pyfftbox;
	PyArrayObject *Pypk;
	PyArrayObject *Pyk;
	PyArrayObject *Pypk2;
	PyArrayObject *Pyk2;

	if(!PyArg_ParseTuple(args, "O!O!O!O!O!", 
		&PyArray_Type, &Pyfftbox, 
		&PyArray_Type, &Pypk, &PyArray_Type, &Pyk,
		&PyArray_Type, &Pypk2, &PyArray_Type, &Pyk2))
		return NULL;
	fft->dim = Pyfftbox->nd;
	//printf("%d \n", fft->dim);
	fft->data = (double*)Pyfftbox->data;
	fft->sizex = Pyfftbox->dimensions[0];
	fft->sizey = Pyfftbox->dimensions[1];
	fft->sizez = Pyfftbox->dimensions[2];
	fft->data = (double*)Pyfftbox->data;
	//printf("%d %d %d\n", fft->sizex, fft->sizey, fft->sizez);
	pk->N = Pypk->dimensions[0];
	pk->val = (double*)Pypk->data;
	pk->k = (double*)Pyk->data;
	pk->Np = Pypk2->dimensions[0];
	pk->Nv = Pypk2->dimensions[1];
	pk->val2 = (double*)Pypk2->data;
	pk->k2 = (double*)Pyk2->data;
	makepk(fft, pk);
	return Py_BuildValue("i", 0);
}

//static PyObject * FFTW(PyObject *self, PyObject *args){
//
//	fftDATA *fftdata = (fftDATA*)malloc(sizeof(fftDATA));
//	PyArrayObject *Pyfftbox;
//
//	if(!PyArg_ParseTuple(args, "O!", &PyArray_Type, &Pyfftbox))
//		return NULL;
//	fftdata->dim = Pyfftbox->nd;
//	fftdata->sizex = Pyfftbox->dimensions[0];
//	fftdata->sizey = Pyfftbox->dimensions[1];
//	fftdata->sizez = Pyfftbox->dimensions[2];
//	fftdata->data = (double*)Pyfftbox->data;
//
//	return Py_BuildValue("i", 0);
//}

static PyMethodDef MakePowerMethods[] = {
//	Define the Method of module
//	{"name used in python", name of the function called in c, something else}
	{"Filling", Filling, METH_VARARGS},
	{"nFilling", nFilling, METH_VARARGS},
	{"Make", Make, METH_VARARGS},
	//{"FFT", FFTW, METH_VARARGS},
	{NULL, NULL}
};

void initMakePower(){
//	m = Py_InitModule("Python Module Name", The Methods of such Module);
// 	The Methods of Module MUST be the same as the name in PyMethodDef
	(void) Py_InitModule("MakePower", MakePowerMethods);
	import_array();// MUST be present for numpy. Called first after above line.
}

//////////////////////////////////////////////////////////////////////////////
