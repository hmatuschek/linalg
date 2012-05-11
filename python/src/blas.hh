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

/** Interface to BLAS DAXPY. */
static PyObject *linalg_blas_daxpy(PyObject *self, PyObject *args);

/** Interface to BLAS DSCAL. */
static PyObject *linalg_blas_dscal(PyObject *self, PyObject *args);

/** Interface to BLAS DGEMV. */
static PyObject *linalg_blas_dgemv(PyObject *self, PyObject *args);

/** Interface to BLAS DTRMV. */
static PyObject *linalg_blas_dtrmv(PyObject *self, PyObject *args);

/** Interface to BLAS DGEMM. */
static PyObject *linalg_blas_dgemm(PyObject *self, PyObject *args);

}

#endif // __LINALG_PYHON_MODULE_HH__
