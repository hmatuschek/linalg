/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */


/**
 * @defgroup python Simple interface to Python/NumPy arrays
 */


#ifndef __LINALG_PYTHON_HH__
#define __LINALG_PYTHON_HH__

#include "matrix.hh"
#include "vector.hh"

extern "C" {
#include <Python.h>
#include <numpy/arrayobject.h>
}


namespace Linalg {

template<>
class NumpyArrayMngr<double> : public DataMngr<double>
{
private:
  /**
   * Holds a reference to the numpy array object.
   */
  PyObject *_array;

public:
  NumpyArrayMngr(PyObject *array)
  throw (PythonError)
    : DataMngr((double *) PyArray_DATA(array), false), _array(array)
  {
    // Check if array is a NumPy array:
    if (0 == PyArray_Check(array)) {
      PythonError err;
      err << "Can not convert to Matrix<double>: Object is not a PyArray instance!";
      throw err;
    }

    // Check if data type is double:
    if (NPY_DOUBLE != PyArray_TYPE(array)) {
      PythonError err;
      err << "Can not convert NumPy array to Matrix<double>: Elements are not of type double.";
      throw err;
    }

    // Check if data is aligned:
    if (0 == PyArray_ISALIGNED(array)) {
      PythonError err;
      err << "Can not convert NumPy array to Matrix<double>: Memory is not aligned!";
      throw err;
    }

    // Check if data is writable:
    if (0 == PyArray_ISWRITEABLE(array))
    {
      PythonError err;
      err << "Can not convert NumPy array to Matrix<double>: Memory is not writable!";
      throw err;
    }

    // Increment reference counter of NumPy array:
    Py_INCREF(_array);
    // Get pointer to data:
    _data = (double *) PyArray_DATA(array);
  }

  /**
   * Destructor. Simply dereferences the array.
   */
  virtual ~PythonError() {
    Py_DECREF(_array);
  }

  /**
   * Retunrs the pointer to the numpy array.
   */
  PyObject *obj() {
    return _array;
  }
};


/**
 * Creates a weak @c Matrix reference from the given NumPy array.
 *
 * The given array must be a 2D aligned NumPy array (PyArray_Type).
 *
 * @ingroup python
 */
inline Matrix<double> doubleMatrixFromNumpyArray(PyObject *array)
{
  NumpyArrayMngr<double> *mngr = new NumpyArrayMngr<double>(array);

  // Check if array is a matrix (2D):
  if (2 != PyArray_NDIM(array)) {
    PythonError err;
    err << "Can not convert NumPy array to Matrix<double>: Array is not of dimension 2!";
    throw err;
  }

  bool is_transposed = false;
  void *data = PyArray_DATA(array);

  size_t rows = PyArray_DIM(array, 0);
  size_t cols = PyArray_DIM(array, 1);

  size_t stride;
  bool   is_rowmajor;
  size_t row_stride = PyArray_STRIDE(array, 0)/NPY_SIZEOF_DOUBLE;
  size_t col_stride = PyArray_STRIDE(array, 1)/NPY_SIZEOF_DOUBLE;

  if (1 == col_stride) {
    is_rowmajor = true;
    stride      = row_stride;
  } else if (1 == row_stride) {
    is_rowmajor = false;
    stride      = col_stride;
  }

  return Matrix<double>::fromData(static_cast<double *>(data), rows, cols, stride, 0,
                                  is_transposed, is_rowmajor);
}


/**
 * Creates a @c Vector view as a weak reference to a 1D NumPy array.
 *
 * @ingroup python
 */
inline Vector<double> doubleVectorFromNumpyArray(PyObject *array)
{
  // Check if array is a NumPy array:
  if (0 == PyArray_Check(array)) {
    PythonError err;
    err << "Can not convert to Vector<double>: Object is not a PyArray instance!";
    throw err;
  }

  // Check if data type is double:
  if (NPY_DOUBLE != PyArray_TYPE(array)) {
    PythonError err;
    err << "Can not convert NumPy array to Vector<double>: Elements are not of type double.";
    throw err;
  }

  // Check if data is aligned:
  if (0 == PyArray_ISALIGNED(array)) {
    PythonError err;
    err << "Can not convert NumPy array to Vector<double>: Memory is not aligned!";
    throw err;
  }

  // Check if data is writable:
  if (0 == PyArray_ISWRITEABLE(array))
  {
    PythonError err;
    err << "Can not convert NumPy array to Vector<double>: Memory is not writable!";
    throw err;
  }

  // Check if array is a matrix (1D):
  if (1 != PyArray_NDIM(array)) {
    PythonError err;
    err << "Can not convert NumPy array to Vector<double>: Array is not of dimension 1!";
    throw err;
  }

  void *data = PyArray_DATA(array);
  size_t dim = PyArray_DIM(array, 0);
  size_t stride = PyArray_STRIDE(array, 0)/NPY_SIZEOF_DOUBLE;

  return Vector<double>(static_cast<double *>(data), dim, 0, stride, false);
}


}


#endif
