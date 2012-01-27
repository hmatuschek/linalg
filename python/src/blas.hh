#ifndef __LINALG_PYHON_MODULE_HH__
#define __LINALG_PYHON_MODULE_HH__

#include <Python.h>

#include "linalg.hh"
#include "python.hh"


extern "C" {

/** Interface to BLAS DDOT. */
static PyObject *linalg_blas_ddot(PyObject *self, PyObject *args);

/** Interface to BLAS DNRM2. */
static PyObject *linalg_blas_dnrm2(PyObject *self, PyObject *args);

/** Debug output for matrix type. */
static PyObject *linalg_blas_print_matrix(PyObject *self, PyObject *args);

}

#endif // __LINALG_PYHON_MODULE_HH__
