#include <Python.h>
#include <arrayobject.h>
#include "cMakePower.h"
//#include "cFFTW.h"

static PyObject * Filling(PyObject *self, PyObject *args){

	FillConf *conf = (FillConf *)malloc(sizeof(FillConf));
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
	fillingf(conf);
	return Py_BuildValue("i", 0);
}
static PyObject * nFilling(PyObject *self, PyObject *args){

	FillConf *conf = (FillConf *)malloc(sizeof(FillConf));
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
