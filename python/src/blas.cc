#include "blas.hh"
#include <iostream>

using namespace Linalg;


extern "C" {


static PyObject *
linalg_blas_ddot(PyObject *self, PyObject *args)
{
  PyObject *py_x,*py_y = 0;
  if (!PyArg_ParseTuple(args, "OO", &py_x, &py_y)) {
    return 0;
  }

  // Get matrix from numpy-array:
  Vector<double> x,y;
  try {
    x = doubleVectorFromNumpyArray(py_x);
    y = doubleVectorFromNumpyArray(py_y);

    return PyFloat_FromDouble(Blas::dot(x, y));
  } catch (Exception &err) {
    PyErr_SetString(PyExc_RuntimeError, err.what());
    return 0;
  }

  return Py_None;
}



static PyObject *
linalg_blas_dnrm2(PyObject *self, PyObject *args)
{
  PyObject *py_x = 0;
  if (!PyArg_ParseTuple(args, "O", &py_x)) {
    return 0;
  }

  // Get matrix from numpy-array:
  Vector<double> x;
  try {
    x = doubleVectorFromNumpyArray(py_x);

    return PyFloat_FromDouble(Blas::nrm2(x));
  } catch (Exception &err) {
    PyErr_SetString(PyExc_RuntimeError, err.what());
    return 0;
  }

  return Py_None;
}



static PyObject *
linalg_blas_dgemv(PyObject *self, PyObject *args)
{
  PyObject *py_A, *py_x, *py_y = 0;
  double alpha, beta;
  if (!PyArg_ParseTuple(args, "dOOdO", &alpha, &py_A, &py_x, &beta, &py_y)) {
    return 0;
  }

  Matrix<double> A;
  Vector<double> x, y;
  try {
    // Crate matrix and vector views:
    A = doubleMatrixFromNumpyArray(py_A);
    x = doubleVectorFromNumpyArray(py_x);
    y = doubleVectorFromNumpyArray(py_y);

    // Perform operation in-place
    Blas::gemv(alpha, A, x, beta, y);

  } catch (Exception &err) {
    // On error...
    PyErr_SetString(PyExc_RuntimeError, err.what());
    return 0;
  }

  return Py_None;
}



/*
 * List all methods:
 */
static PyMethodDef LinalgBlasMethods[] = {
  {"ddot", linalg_blas_ddot, METH_VARARGS, "BLAS level-1 DDOT() wrapper."},
  {"dnrm2", linalg_blas_dnrm2, METH_VARARGS, "BLAS level-1 DNRM2() wrapper."},
  {"dgemv", linalg_blas_dgemv, METH_VARARGS, "BLAS level-2 DGEMV() wrapper."},
  {NULL, NULL, 0, NULL}
};



/*
 * Initialize module
 */
PyMODINIT_FUNC
init_blas(void)
{
  (void) Py_InitModule("blas", LinalgBlasMethods);

  import_array();
}

}
