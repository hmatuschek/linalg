/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_MATRIX_HH__
#define __LINALG_MATRIX_HH__

#include "memory.hh"
#include "vector.hh"


namespace Linalg {


/**
 * Defines a matrix.
 *
 * @todo Implement "empty" matrix constructor.
 *
 * @ingroup linalg
 */
template <class Scalar>
class Matrix : public DataPtr<Scalar>
{
public:
  typedef std::auto_ptr< Matrix<Scalar> > unowned;

  /**
   * Represents a weak reference to the values of a matrix. This can be used to copy values of
   * one matrix into an other.
   */
  class ValueRef
  {
  protected:
    Matrix<Scalar> _matrix;

  public:
    ValueRef(const Matrix &matrix)
      : _matrix(matrix)
    {
      // Pass...
    }

    ValueRef(const ValueRef &other)
      : _matrix(other._matrix)
    {
      // Pass...
    }

    const ValueRef &operator= (const Scalar &val)
    {
      for (size_t i=0; i<this->_matrix.rows(); i++) {
        for (size_t j=0; j<this->_matrix.cols(); j++) {
          this->_matrix(i,j) = val;
        }
      }
    }

    const ValueRef &operator= (const Matrix &other)
    {
      // Assert equal shapes:
      LINALG_SHAPE_ASSERT(this->_matrix.rows() == other.rows());
      LINALG_SHAPE_ASSERT(this->_matrix.cols() == other.cols());

      for (size_t i=0; i<this->_matrix.rows(); i++) {
        for (size_t j=0; j<this->_matrix.cols(); j++) {
          this->_matrix(i,j) = other(i,j);
        }
      }
    }
  };


protected:
  /**
   * Holds the number of rows of the matrix.
   */
  size_t _rows;

  /**
   * Holds the number of columns of the matrix.
   */
  size_t _cols;

  /**
   * Number of rows in the outer most matrix if matrix is in column-major form and number of
   * columns if matrix is in row-major form.
   */
  size_t _stride;

  /**
   * Holds the offset of the first element in the matrix, starting from the first element of the
   * original memory element.
   */
  size_t _offset;

  /**
   * If true, the matrix is assumed to be transposed.
   */
  bool _transposed;

  /**
   * If true, the matrix is stored in row-major (C) order, if false then it is stored in
   * column major (Fortran) order.
   */
  bool _is_rowmajor;


protected:
  /**
   * Full constructor. Shall not be used from the outside.
   */
  Matrix(Scalar *data, size_t rows, size_t cols,
         size_t stride, size_t offset,
         bool transposed, bool rowmajor, bool owns_data)
    : DataPtr<Scalar>(data, owns_data), _rows(rows), _cols(cols),
      _stride(stride), _offset(offset),
      _transposed(transposed), _is_rowmajor(rowmajor)
  {
    // Pass...
  }

  /**
   * Calculates the index (in flattened data array) for the given row/col index.
   *
   * @param i Specifies the row index (zero-based).
   * @param j Specifies the column index (zero-based).
   */
  inline size_t _getIndex(size_t i, size_t j) const
  {
    if (this->_transposed)
      std::swap(i,j);

    if (this->_is_rowmajor)
      return this->_offset + this->_stride*i + j;
    else
      return this->_offset + this->_stride*j + i;
  }


public:
  /**
   * Constructs an empty matrix.
   */
  Matrix()
    : DataPtr<Scalar>(), _rows(0), _cols(0), _stride(0), _offset(0),
      _transposed(false), _is_rowmajor(true)
  {
    // Pass...
  }


  /**
   * Copy constructor with ownership transfer.
   */
  explicit Matrix(Matrix<Scalar> &other, bool take_ownership)
    : DataPtr<Scalar>(other, take_ownership), _rows(other._rows), _cols(other._cols),
      _stride(other._stride), _offset(other._offset),
      _transposed(other._transposed), _is_rowmajor(other._is_rowmajor)
  {
    // Pass...
  }


  /**
   * Takes the ownership of a matrix that is unowned.
   */
  explicit Matrix(unowned other)
    : DataPtr<Scalar>(*other, true), _rows(other->_rows), _cols(other->_cols),
      _stride(other->_stride), _offset(other->_offset),
      _transposed(other->_transposed), _is_rowmajor(other->_is_rowmajor)
  {
    // explicity take ownership
    other.release();
  }

  /**
   * Copy constructor with weak reference.
   */
  Matrix(const Matrix<Scalar> &other)
    : DataPtr<Scalar>(other), _rows(other._rows), _cols(other._cols),
      _stride(other._stride), _offset(other._offset),
      _transposed(other._transposed), _is_rowmajor(other._is_rowmajor)
  {
    // Pass...
  }

  /**
   * Constructs a new matrix (allocates new memory) with given rows and columns.
   */
  Matrix(size_t rows, size_t cols, bool rowmajor=true)
    : DataPtr<Scalar>(rows*cols), _rows(rows), _cols(cols),
      _stride(cols), _offset(0), _transposed(false), _is_rowmajor(rowmajor)
  {
    if (! this->_is_rowmajor)
      this->_stride = this->_rows;
  }

  /**
   * Assignment operator for an unowned matrix.
   *
   * This matrix takes the ownership of the data of the assigned matrix.
   */
  inline const Matrix<Scalar> &operator= (Matrix<Scalar>::unowned other)
  {
    // Free data if owning it:
    if(this->_owns_data && this->_data != other->_data)
    {
      delete this->_data;
    }

    // Take ownership of data:
    this->_data = other->_data; other->releaseData();
    this->_owns_data = true;

    // Copy meta-data
    this->_rows = other->_rows;
    this->_cols = other->_cols;
    this->_stride = other->_stride;
    this->_transposed = other->_transposed;
    this->_is_rowmajor = other->_is_rowmajor;

    return *this;
  }


  /**
   * Assignment (weak) operator.
   */
  inline const Matrix<Scalar> &operator= (const Matrix<Scalar> &other)
  {
    if (this->_owns_data && other._data != this->_data) {
      delete this->_data;
    }

    this->_data = other._data;
    this->_owns_data = false;

    this->_rows = other._rows;
    this->_cols = other._cols;
    this->_stride = other._stride;
    this->_transposed = other._transposed;
    this->_is_rowmajor = other._is_rowmajor;

    return *this;
  }


  /**
   * Explicit copy of the matrix.
   *
   * This method returns a new (unowned) matrix.
   */
  inline unowned copy()
  {
    // Construct dense row-major matrix:
    Matrix ret(this->rows(), this->cols(), this->isRowOrder());

    // copy values:
    for (size_t i=0; i<this->rows(); i++)
    {
      for (size_t j=0; j<this->cols(); j++)
      {
        ret(i,j) = (*this)(i,j);
      }
    }

    // Done...
    return ret.takeOwnership();
  }


  inline ValueRef vals()
  {
    return ValueRef(*this);
  }


  /**
   * Returns the number of rows of the matrix.
   *
   * @note This method does not return the number of rows in the memory representation!
   */
  inline size_t rows() const
  {
    if (this->_transposed)
      return this->_cols;
    return this->_rows;
  }


  /**
   * Returns the number of columns.
   *
   * @note This method does not return the number of columns in the memory representation!
   */
  inline size_t cols() const
  {
    if (this->_transposed)
      return this->_rows;
    return this->_cols;
  }


  /**
   * Returns the number of rows if matrix in in column-major (Fortran) form otherwise
   * return the number columns of the outer-most matrix.
   */
  inline size_t stride() const
  {
    return this->_stride;
  }


  /**
   * Retuns true, if the matrix is transposed.
   */
  inline bool isTransposed() const {
    return this->_transposed;
  }


  /**
   * Returns true, if the storage-order of the matrix is row-major (C order).
   */
  inline bool isRowOrder() const {
    return this->_is_rowmajor;
  }


  /**
   * Returns a matrix in column-major order.
   */
  inline Matrix<Scalar> colMajor() const
  {
    // If this matrix is in row-major:
    if (this->_is_rowmajor)
    {
      return Matrix<Scalar>(this->_data, this->_cols, this->_rows,
                            this->_stride, this->_offset,
                            ! this->isTransposed(), false,
                            false);
    }

    // If matrix is in column-major:
    return *this;
  }


  /**
   * Returns a matrix in column-major order.
   */
  inline Matrix<Scalar> rowMajor() const
  {
    // If this matrix is in col-major:
    if (! this->_is_rowmajor)
    {
      return Matrix<Scalar>(this->_data,
                            this->_cols, this->_rows, this->_stride, this->_offset,
                            ! this->isTransposed(), true,
                            false);
    }

    // If matrix is in column-major:
    return *this;
  }


  /**
   * Returns the pointer to the first element of the array.
   */
  inline Scalar* operator* ()
  {
    return this->_data + this->_offset;
  }


  /**
   * Returns a const pointer to the first element of the array.
   */
  inline const Scalar* operator* () const
  {
    return this->_data + this->_offset;
  }


  /**
   * Returns a reference to the element of the i-th row and j-th column in the matrix.
   */
  inline Scalar &operator() (size_t i, size_t j)
  {
    return this->_data[this->_getIndex(i, j)];
  }


  /**
   * Returns a const reference to the element of the i-th row and j-th column in the matrix.
   */
  const Scalar &operator() (size_t i, size_t j) const
  {
    return this->_data[this->_getIndex(i,j)];
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
    if (this->_transposed)
    {
      if (i >= this->_cols)
      {
        Linalg::IndexError err;
        err << "Row-index " << i << " out of range 0-" << this->_cols-1;
        throw err;
      }

      size_t stride = this->stride();
      if (! this->_is_rowmajor)
        stride = 1;

      return Vector<Scalar>(this->_data, this->_rows, this->_getIndex(i,0), stride, false);
    }

    if (i >= this->_rows)
    {
      Linalg::IndexError err;
      err << "Row-index " << i << " out of range 0-" << this->_rows-1;
      throw err;
    }

    size_t stride = 1;
    if (! this->_is_rowmajor)
      stride = this->_stride;

    // Construct vector from matrix
    return Vector<Scalar>(this->_data, this->_cols, this->_getIndex(i, 0), stride, false);
  }


  /**
   * Returns the j-th column as a vector.
   */
  inline Vector<Scalar> col(size_t j)
  {
    if (this->_transposed)
    {
      if (j >= this->_rows)
      {
        Linalg::IndexError err;
        err << "Row-index " << j << " out of range 0-" << this->_rows-1;
        throw err;
      }

      size_t stride = 1;
      if (! this->_is_rowmajor)
        stride = this->_stride;

      // Construct vector from matrix
      return Vector<Scalar>(this->_data, this->_cols, this->_getIndex(0, j), stride, false);
    }

    if (j >= this->_cols)
    {
      Linalg::IndexError err;
      err << "Row-index " << j << " out of range 0-" << this->_cols-1;
      throw err;
    }

    size_t stride = this->_stride;
    if (! this->_is_rowmajor)
      stride = 1;

    return Vector<Scalar>(this->_data, this->_rows, this->_getIndex(0, j), stride, false);
  }


  /**
   * Returns a sub-matrix (block) of the matrix.
   */
  inline Matrix<Scalar> sub(size_t i, size_t j, size_t rows, size_t cols)
  {
    // Test if row-indices match
    if (i >= this->_rows || i+rows > this->_rows)
    {
      Linalg::IndexError err;
      err << "Row indices " << i << "-" << i+rows-1 << " not in range 0-" << this->_rows-1;
      throw err;
    }

    // Test if column-indices match:
    if (j >= this->_cols || j+cols > this->_cols)
    {
      Linalg::IndexError err;
      err << "Column indices " << j << "-" << j+cols-1 << " not in range 0-" << this->_cols-1;
      throw err;
    }

    // Assemble new Matrix<T> view:
    return Matrix<Scalar>(this->_data,
                          rows, cols, this->_stride, this->_getIndex(i,j),
                          this->_transposed, this->_is_rowmajor,
                          false);
  }


  /**
   * Returns a transposed view of the matrix.
   */
  inline Matrix<Scalar> t() const
  {
    return Matrix<Scalar>(this->_data,
                          this->_rows, this->_cols, this->_stride, this->_offset,
                          !this->_transposed, this->_is_rowmajor, false);
  }


  inline unowned takeOwnership()
  {
    return unowned(new Matrix<Scalar>(*this, true));
  }


public:
  /**
   * Constructs a @c Matrix from the given data.
   */
  static Matrix<Scalar> fromData(Scalar *data, size_t rows, size_t cols,
                                 size_t stride=0, size_t offset=0,
                                 bool transposed=false, bool row_major=true)
  {
    size_t _stride = cols;
    if (0 != stride) {
      _stride = stride;
    } else if (! row_major) {
      _stride = rows;
    }

    return Matrix<Scalar>(data, rows, cols, _stride, offset, transposed, row_major, false);
  }


  static unowned
  fromUnownedData(Scalar *data, size_t rows, size_t cols, size_t stride=0, size_t offset=0,
                  bool transposed=false, bool row_major=true)
  {
    size_t _stride = cols;
    if (0 != stride) {
      _stride = stride;
    } else if (! row_major) {
      _stride = rows;
    }

    return Matrix<Scalar>(data, rows, cols, _stride, offset, transposed, row_major, false).takeOwnership();
  }


  /**
   * Returns a matrix initialized with all values = 0.
   */
  static unowned zeros(size_t rows, size_t cols)
  {
    Matrix<Scalar> m(rows, cols);

    for (size_t i=0; i<rows; i++)
    {
      for (size_t j=0; j<cols; j++)
      {
        m(i,j) = 0;
      }
    }

    return m.takeOwnership();
  }


  /**
   * Returns a unit matrix.
   */
  static unowned unit(size_t rows, size_t cols)
  {
    Matrix<Scalar> U = Matrix<Scalar>::zeros(rows, cols);

    for (size_t i=0; i<std::min(rows, cols); i++)
    {
      U(i,i) = 1;
    }

    return U.takeOwnership();
  }


  /**
   * Returns a squere-unit matrix.
   */
  static unowned unit(size_t N)
  {
    return Matrix<Scalar>::unit(N,N);
  }


  /**
   * Constructs a new matrix initialized with [0,1] uniform distributed random-numbers.
   */
  static unowned rand(size_t rows, size_t cols)
  {
    Matrix<Scalar> m(rows, cols);

    for (size_t i=0; i<rows; i++)
    {
      for (size_t j=0; j<cols; j++)
      {
        m(i,j) = Scalar(std::rand())/RAND_MAX;
      }
    }

    return m.takeOwnership();
  }
};


}

#endif
