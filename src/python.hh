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

/**
 * Template declaration for the python numpy-array data manager.
 */
template <class Scalar> class NumpyArrayMngr;

/**
 * Implementation of the python numpy-array data manager (doubles).
 *
 * This class manages a reference to a numpy-array object and extracts the pointer
 * to the data. On construction, the reference counter of the numpy array is increased and
 * decreased on destruction of the manager.
 *
 * @ingroup python
 */
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
    : DataMngr<double>((double *) PyArray_DATA(array), false), _array(array)
  {
    // Increment reference counter of NumPy array:
    Py_INCREF(_array);

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
  }


  /**
   * Destructor. Simply dereferences the array.
   */
  virtual ~NumpyArrayMngr() {
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
  DataPtr<double> data(new NumpyArrayMngr<double>(array));

  // Check if array is a matrix (2D):
  if (2 != PyArray_NDIM(array)) {
    PythonError err;
    err << "Can not convert NumPy array to Matrix<double>: Array is not of dimension 2!";
    throw err;
  }

  size_t rows = PyArray_DIM(array, 0);
  size_t cols = PyArray_DIM(array, 1);
  size_t row_stride = PyArray_STRIDE(array, 0)/NPY_SIZEOF_DOUBLE;
  size_t col_stride = PyArray_STRIDE(array, 1)/NPY_SIZEOF_DOUBLE;

  return Matrix<double>(data, 0, rows, cols, row_stride, col_stride);
}


/**
 * Creates a @c Linalg::Vector view from reference to a 1D NumPy array.
 *
 * @ingroup python
 */
inline Vector<double> doubleVectorFromNumpyArray(PyObject *array)
{
  // Create data-reference from array object:
  DataPtr<double> data(new NumpyArrayMngr<double>(array));

  // Check if array is a vector (1D):
  if (1 != PyArray_NDIM(array)) {
    PythonError err;
    err << "Can not convert NumPy array to Vector<double>: Array is not of dimension 1!";
    throw err;
  }

  size_t dim = PyArray_DIM(array, 0);
  size_t stride = PyArray_STRIDE(array, 0)/NPY_SIZEOF_DOUBLE;
  return Vector<double>(data, dim, 0, stride);
}


}


#endif
