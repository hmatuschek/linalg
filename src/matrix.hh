#ifndef __LINALG_MATRIX_HH__
#define __LINALG_MATRIX_HH__


namespace Linalg {


/**
 * Defines a view to a matrix.
 *
 * @ingroup linalg
 */
template <class T>
class Matrix
{
protected:
  /**
   * Holds a managed reference to the data array, the matrix lives in.
   */
  SmartPtr< ArrayBase<T> > _data;

  /**
   * Holds the number of rows of the matrix.
   */
  size_t _rows;

  /**
   * Holds the number of columns of the matrix.
   */
  size_t _cols;

  /**
   * Holds the leading dimension of the matrix.
   *
   * This is used to address a matrix-in-matrix. Usually, the LD (leading dimension) is the number
   * of rows of the most-outer matrix.
   */
  size_t _leading_dimension;

  size_t _offset;

  bool _transposed;


protected:
  /**
   * Full constructor. Shall not be used from the outside.
   */
  Matrix(SmartPtr< ArrayBase<T> > data, size_t rows, size_t cols, size_t ld, size_t offset, bool transposed)
    : _data(data), _rows(rows), _cols(cols), _leading_dimension(ld), _offset(offset), _transposed()
  {
    // Pass...
  }

  /**
   * Calculates the index (in flattened data array) for the given row/col index.
   */
  inline size_t _getIndex(size_t i, size_t j) const
  {
    return this->_offset + this->_leading_dimension*i + j;
  }


public:
  /**
   * Copy constructor, does not copy the memory of the matrix.
   */
  Matrix(const Matrix<T> &other)
    : _data(other._data), _rows(other._rows), _cols(other._cols),
      _leading_dimension(other._leading_dimension), _offset(other._offset)
  {
    // Pass...
  }


  /**
   * Constructs a new matrix (allocates new memory) with given rows and columns.
   */
  Matrix(size_t rows, size_t cols)
    : _data(0), _rows(rows), _cols(cols), _leading_dimension(rows), _offset(0)
  {
    // allocate some data
    this->_data = SmartPtr< ArrayBase<T> >(new ArrayBase<T>(rows*cols));
  }


  /**
   * True copy, also copies the memory.
   */
  Matrix<T> copy()
  {
    return Matrix<T>(SmartPtr< ArrayBase<T> >(new ArrayBase<T>(*(this->_data.get_ptr()))),
                     this->_rows, this->_cols, this->_leading_dimension, this->_offset);
  }


  inline size_t rows() const
  {
    return this->_rows;
  }


  inline size_t cols() const
  {
    return this->_cols;
  }


  inline size_t leading_dimension() const
  {
    return this->_leading_dimension;
  }


  inline T* operator* ()
  {
    return this->_data->getPtr() + sizeof(T)*this->_offset;
  }


  const inline T* operator* () const
  {
    return this->_data->getPtr() + sizeof(T)*this->_offset;
  }


  T &operator() (size_t i, size_t j)
  {
    return this->_data->getPtr()[this->_getIndex(i, j)];
  }


  const T &operator() (size_t i, size_t j) const
  {
    return this->_data->getPtr()[this->_getIndex(i,j)];
  }


  Vector<T> row(size_t i)
  {
    if (i >= this->_rows)
    {
      Linalg::IndexError err;
      err << "Row-index " << i << " out of range 0-" << this->_rows-1;
    }

    // Construct vector from matrix
    Vector<T>(this->_data, this->_rows, this->_getIndex(i, 0), this->_leading_dimension);
  }


  Vector<T> col(size_t j)
  {
    if (j >= this->_cols)
    {
      Linalg::IndexError err;
      err << "Column-index " << j << " out of range 0-"<< this->_cols-1;
    }
  }


  Matrix<T> sub(size_t i, size_t j, size_t rows, size_t cols)
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
    return Matrix<T>(this->_data, rows, cols, this->_leading_dimension, this->_getIndex(i,j));
  }


  void set(const Matrix<T> &other)
  {
    // Check dims:
    if (this->rows() != other.rows() || this->cols() != other.cols())
    {
      throw Linalg::IndexError("Shape mismatch!");
    }

    for (size_t i=0; i<this->rows(); i++)
    {
      for (size_t j=0; j<this->cols(); j++)
      {
        (*this)(i,j) = other(i,j);
      }
    }
  }


  inline void set(size_t i, size_t j, size_t rows, size_t cols, const Matrix<T> other)
  {
    this->sub(i,j, rows,cols).set(other);
  }

  
  inline Matrix<T> T() const
  {
    return Matrix<T>(this->_data, this->_rows, this->_columns, this->_leading_dimension, this->_offset, !this->_transposed);
  }



public:
  static Matrix<T> zeros(size_t rows, size_t cols)
  {
    Matrix<T> m(rows, cols);

    for (size_t i=0; i<rows; i++)
    {
      for (size_t j=0; j<cols; j++)
      {
        m(i,j) = 0;
      }
    }

    return m;
  }


  static Matrix<T> unit(size_t rows, size_t cols)
  {
    Matrix<T> U = Matrix<T>::zeros(rows, cols);

    for (size_t i=0; i<std::min(rows, cols); i++)
    {
      U(i,i) = 1;
    }

    return U;
  }


  static Matrix<T> unit(size_t N)
  {
    return Matrix<T>::unit(N,N);
  }


  static Matrix<T> rand(size_t rows, size_t cols)
  {
    Matrix<T> m(rows, cols);

    for (size_t i=0; i<rows; i++)
    {
      for (size_t j=0; j<cols; j++)
      {
        m(i,j) = T(std::rand())/RAND_MAX;
      }
    }

    return m;
  }
};


}

#endif
