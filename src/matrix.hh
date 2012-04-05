/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_MATRIX_HH__
#define __LINALG_MATRIX_HH__

#include "array.hh"
#include "vector.hh"
#include "blas/utils.hh"
#include <iostream>


namespace Linalg {


/**
 * Defines a matrix.
 *
 * @ingroup matrix
 */
template <class Scalar>
class Matrix
    : public Array<Scalar>
{
protected:
  /**
   * Calculates the index (in flattened data array) for the given row/col index.
   *
   * @param i Specifies the row index (zero-based).
   * @param j Specifies the column index (zero-based).
   */
  inline size_t _getIndex(size_t i, size_t j) const
  {
    return this->_offset + i*this->_strides[0] + j*this->_strides[1];
  }


public:
  /**
   * Constructor for an empty matrix.
   */
  Matrix(size_t rows, size_t cols, bool rowmajor=true)
    : Array<Scalar>()
  {
    this->_data = new SmartPtr<Scalar>(new Scalar[rows*cols], true);
    this->_offset = 0;
    this->_shape.resize(2); this->_shape[0] = rows; this->_shape[1] = cols;
    this->_strides.resize(2);
    if (rowmajor) {
      this->_strides[0] = cols;
      this->_strides[1] = 1;
    } else {
      this->_strides[0] = 1;
      this->_strides[1] = rows;
    }
  }


  /**
   * Constructs an empty matrix.
   */
  Matrix()
    : Array<Scalar>()
  {
    // Pass...
  }


  /**
   * Copy constructor with ownership transfer.
   */
  Matrix(const Matrix<Scalar> &other)
    : Array<Scalar>(other._data, other._offset, other._shape, other._strides)
  {
    // Pass...
  }

  /**
   * Copy constructor with ownership transfer.
   */
  Matrix(const Array<Scalar> &other)
    : Array<Scalar>(other.data(), other.offset(), other.shape(), other.strides())
  {
    LINALG_SHAPE_ASSERT(2 == other.ndim());
  }


  /**
   * Explicit copy of the matrix.
   *
   * This method returns a new (unowned) matrix.
   */
  inline Matrix<Scalar> copy()
  {
    return Matrix<Scalar>(Array<Scalar>::copy());
  }

  void swap(Matrix<Scalar> &other)
  {
    Matrix<Scalar> tmp = other;
    other = *this;
    Array<Scalar>::operator =(tmp);
  }

  Matrix<Scalar> t() const {
    return Array<Scalar>::t();
  }


  /**
   * Returns the number of rows of the matrix.
   *
   * @note This method does not return the number of rows in the memory representation!
   */
  inline size_t rows() const
  {
    return this->_shape[0];
  }


  /**
   * Returns the number of columns.
   *
   * @note This method does not return the number of columns in the memory representation!
   */
  inline size_t cols() const
  {
    return this->_shape[1];
  }


  /**
   * Returns a reference to the element of the i-th row and j-th column in the matrix.
   */
  inline Scalar &operator() (size_t i, size_t j)
  {
    return this->_data->ptr()[this->_getIndex(i, j)];
  }


  /**
   * Returns a const reference to the element of the i-th row and j-th column in the matrix.
   */
  const Scalar &operator() (size_t i, size_t j) const
  {
    return this->_data->ptr()[this->_getIndex(i,j)];
  }


  /**
   * Returns the (0,0) element of the matrix if the matrix is 1,1 dimensional.
   */
  inline Scalar asScalar() const {
    if (1 != this->_rows || 1 != this->_cols)
    {
      IndexError err; err << "Can not convert matrix to scalar, matrix has more than one element.";
      throw err;
    }

    return (*this)(0,0);
  }


  /**
   * Returns the i-th row as a vector.
   */
  inline Vector<Scalar> row(size_t i)
  {
    return Vector<Scalar>(this->_data, this->_getIndex(i,0), rows(), Array<Scalar>::strides(1));
  }


  /**
   * Returns the j-th column as a vector.
   */
  inline Vector<Scalar> col(size_t j)
  {
    return Vector<Scalar>(this->_data, this->_getIndex(0,j), cols(), Array<Scalar>::strides(0));
  }


  /**
   * Returns a sub-matrix (block) of the matrix.
   */
  inline Matrix<Scalar> sub(size_t i, size_t j, size_t nrow, size_t ncol)
  {
    // Test if row-indices match
    if (i >= rows() || i+nrow > rows())
    {
      Linalg::IndexError err;
      err << "Row indices " << i << "-" << i+nrow-1 << " not in range 0-" << rows()-1;
      throw err;
    }

    // Test if column-indices match:
    if (j >= cols() || j+ncol > cols())
    {
      Linalg::IndexError err;
      err << "Column indices " << j << "-" << j+ncol-1 << " not in range 0-" << cols()-1;
      throw err;
    }

    std::vector<size_t> shape(2);
    shape[0] = nrow; shape[1] = ncol;
    return Matrix<Scalar>(this->_data, this->getIndex(i,j), shape, this->strides());
  }


public:
  /**
   * Constructs a @c Matrix from the given data.
   */
  static Matrix<Scalar> fromData(Scalar *data,
                                 size_t rows, size_t cols, size_t rstride, size_t cstride,
                                 size_t offset=0)
  {
    DataPtr<Scalar> *data_ptr = new SmartPtr<Scalar>(data, false);

    std::vector<size_t> dims(2); dims[0]=rows; dims[1]=cols;
    std::vector<size_t> strd(2); strd[0]=rstride; strd[1]=cstride;
    Matrix<Scalar> ret(Array<Scalar>(data_ptr, offset, dims, strd));

    data_ptr->unref();
    return ret;
  }


  /**
   * Constructs a matrix from the given data and takes the ownership of the data.
   */
  static Matrix<Scalar>
  fromUnownedData(Scalar *data, size_t rows, size_t cols, size_t rstride, size_t cstride,
                  size_t offset=0)
  {
    std::vector<size_t> dims(2); dims[0]=rows; dims[1]=cols;
    std::vector<size_t> strd(2); strd[0]=rstride; strd[1]=cstride;
    DataPtr<Scalar> *data_ptr(data, true);
    Matrix<Scalar> ret(Array<Scalar>(data_ptr, offset, dims, strd));
    data_ptr->unref();
    return ret;
  }


  /**
   * Constructs a matrix with uninitilized values.
   */
  static Matrix<Scalar> empty(size_t rows, size_t cols, bool rowmajor=true)
  {
    DataPtr<Scalar> *data = new SmartPtr<Scalar>(new Scalar[rows*cols], true);
    std::vector<size_t> shape(2); shape[0] = rows; shape[1] = cols;
    std::vector<size_t> strides(2);

    if (rowmajor) {
      strides[0] = cols; strides[1] = 1;
    } else {
      strides[0] = 1; strides[1] = rows;
    }

    Matrix<Scalar> ret(Array<Scalar>(data, 0, shape, strides)); data->unref();

    return ret;
  }


  /**
   * Returns a matrix initialized with all values = 0.
   */
  static Matrix<Scalar> zeros(size_t rows, size_t cols)
  {
    Matrix<Scalar> m(rows, cols);

    for (size_t i=0; i<rows; i++)
    {
      for (size_t j=0; j<cols; j++)
      {
        m(i,j) = 0;
      }
    }

    return m;
  }


  /**
   * Returns a unit matrix.
   */
  static Matrix<Scalar> unit(size_t rows, size_t cols)
  {
    Matrix<Scalar> U = Matrix<Scalar>::zeros(rows, cols);

    for (size_t i=0; i<std::min(rows, cols); i++)
    {
      U(i,i) = 1;
    }

    return U;
  }


  /**
   * Returns a squere-unit matrix.
   */
  static Matrix<Scalar> unit(size_t N)
  {
    return Matrix<Scalar>::unit(N,N);
  }


  /**
   * Constructs a new matrix initialized with [0,1] uniform distributed random-numbers.
   */
  static Matrix<Scalar> rand(size_t rows, size_t cols)
  {
    Matrix<Scalar> m(rows, cols);

    for (size_t i=0; i<rows; i++)
    {
      for (size_t j=0; j<cols; j++)
      {
        m(i,j) = Scalar(std::rand())/RAND_MAX;
      }
    }

    return m;
  }
};

}

#endif
