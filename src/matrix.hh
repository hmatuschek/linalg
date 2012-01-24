#ifndef __LINALG_MATRIX_HH__
#define __LINALG_MATRIX_HH__

#include "arraybase.hh"
#include "vector.hh"


namespace Linalg {


/**
 * Defines a matrix.
 *
 * @ingroup linalg
 */
template <class Scalar>
class Matrix
{
protected:
  /**
   * Holds a managed reference to the data array, the matrix lives in.
   */
  ArrayBase<Scalar> _data;

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
  Matrix(ArrayBase<Scalar> &data, size_t rows, size_t cols,
         size_t stride, size_t offset,
         bool transposed, bool rowmajor)
    : _data(data), _rows(rows), _cols(cols),
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
   * Copy constructor, does not copy the memory of the matrix.
   */
  Matrix(Matrix<Scalar> &other)
    : _data(other._data), _rows(other._rows), _cols(other._cols),
      _stride(other._stride), _offset(other._offset),
      _transposed(other._transposed), _is_rowmajor(other._is_rowmajor)
  {
    // Pass...
  }

  Matrix(const Matrix<Scalar> &other)
    : _data(other._data.weak()), _rows(other._rows), _cols(other._cols),
      _stride(other._stride), _offset(other._offset),
      _transposed(other._transposed), _is_rowmajor(other._is_rowmajor)
  {
    // Pass...
  }

  /**
   * Constructs a new matrix (allocates new memory) with given rows and columns.
   */
  Matrix(size_t rows, size_t cols, bool rowmajor=true)
    : _data(rows*cols), _rows(rows), _cols(cols),
      _stride(cols), _offset(0), _transposed(false), _is_rowmajor(rowmajor)
  {
    if (! this->_is_rowmajor)
      this->_stride = this->_rows;
  }

  /**
   * Assignment operator.
   */
  const Matrix<Scalar> &operator= (Matrix<Scalar> &other)
  {
    this->_data = other._data;
    this->_rows = other._rows;
    this->_cols = other._cols;
    this->_stride = other._stride;
    this->_transposed = other._transposed;
    this->_is_rowmajor = other._is_rowmajor;

    return *this;
  }


  /**
   * Assignment (weak) operator.
   */
  const Matrix<Scalar> &operator= (const Matrix<Scalar> &other)
  {
    this->_data = other._data.weak();
    this->_rows = other._rows;
    this->_cols = other._cols;
    this->_stride = other._stride;
    this->_transposed = other._transposed;
    this->_is_rowmajor = other._is_rowmajor;

    return *this;
  }


  /**
   * Explicit copy.
   */
  inline Matrix<Scalar> copy()
  {
    // Construct dense row-major matrix:
    Matrix ret(this->rows(), this->cols());
    // copy values:
    for (size_t i=0; i<this->rows(); i++)
    {
      for (size_t j=0; j<this->cols(); j++)
      {
        ret(i,j) = (*this)(i,j);
      }
    }
    // Done...
    return ret;
  }

  /**
   * Returns a weak reference to this array.
   */
  inline Matrix<Scalar> weak() const
  {
    ArrayBase<Scalar> data(this->_data.weak());
    return Matrix<Scalar>(data, this->_rows, this->_cols,
                          this->_stride, this->_offset,
                          this->_transposed, this->_is_rowmajor);
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
   * Returns a matrix in column-major order.
   */
  inline Matrix<Scalar> colMajor() const
  {
    // If this matrix is in row-major:
    if (this->_is_rowmajor)
    {
      ArrayBase<Scalar> data(this->_data.weak());
      return Matrix<Scalar>(data, this->_cols, this->_rows,
                            this->_stride, this->_offset,
                            ! this->isTransposed(), false);
    }

    // If matrix is in column-major:
    return this->weak();
  }

  /**
   * Returns a matrix in column-major order.
   */
  inline Matrix<Scalar> rowMajor() const
  {
    // If this matrix is in col-major:
    if (! this->_is_rowmajor)
    {
      return Matrix<Scalar>(this->_data.weak(), this->_cols, this->_rows,
                            this->_stride, this->_offset,
                            ! this->isTransposed(), true);
    }

    // If matrix is in column-major:
    return this->weak();
  }

  /**
   * Returns the pointer to the first element of the array.
   */
  inline Scalar* operator* ()
  {
    return *(this->_data) + this->_offset;
  }

  /**
   * Returns a const pointer to the first element of the array.
   */
  inline const Scalar* operator* () const
  {
    return *(this->_data) + this->_offset;
  }

  /**
   * Returns a reference to the element of the i-th row and j-th column in the matrix.
   */
  Scalar &operator() (size_t i, size_t j)
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

      return Vector<Scalar>(this->_data.weak(), this->_rows, this->_getIndex(i,0), stride);
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
    return Vector<Scalar>(this->_data.weak(), this->_cols, this->_getIndex(i, 0), stride);
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
      return Vector<Scalar>(this->_data.weak(), this->_cols, this->_getIndex(0, j), stride);
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

    return Vector<Scalar>(this->_data.weak(), this->_rows, this->_getIndex(0, j), stride);
  }

  /**
   * Returns a sub-matrix (block) of the matrix.
   */
  Matrix<Scalar> sub(size_t i, size_t j, size_t rows, size_t cols)
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
    return Matrix<Scalar>(this->_data.weak(), rows, cols, this->_stride, this->_getIndex(i,j));
  }

  /**
   * Returns a transposed view of the matrix.
   */
  inline Matrix<Scalar> t() const
  {
    ArrayBase<Scalar> data(this->_data.weak());
    return Matrix<Scalar>(data, this->_rows,
                          this->_cols, this->_stride,
                          this->_offset, !this->_transposed, this->_is_rowmajor);
  }


public:
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
