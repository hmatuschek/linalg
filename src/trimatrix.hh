#ifndef __LINALG_TRIMATRIX_HH__
#define __LINALG_TRIMATRIX_HH__

#include "matrix.hh"

namespace Linalg {

/**
 * Implements a view on @c Matrix as a unit or non-unit upper- or lower-triangular matrix.
 */
template <class Scalar>
class TriMatrix : public Matrix<Scalar>
{
protected:
  /**
   * Specifies if the triangular matrix is storred in the upper or lower half.
   */
  bool _is_upper;

  /**
   * Specifies if the matrix has a unit diagonal.
   */
  bool _is_unit_triangular;


public:
  /**
   * Constructs a triangular matrix-view from the given marix.
   */
  TriMatrix(Matrix<Scalar> &matrix, bool upper, bool unit)
    : Matrix<Scalar>(matrix.weak()), _is_upper(upper), _is_unit_triangular(unit)
  {
    if (matrix.transposed())
      this->_is_upper = ! this->_is_upper;
  }


  /**
   * Copy constructor.
   */
  TriMatrix(TriMatrix<Scalar> &other)
    : Matrix<Scalar>(other.weak()),
      _is_upper(other._is_upper), _is_unit_triangular(other._is_unit_triangular)
  {
    // Pass...
  }


  /**
   * Assignment operation.
   */
  TriMatrix<Scalar> &operator= (TriMatrix<Scalar> &other)
  {
    // Assign Matrix:
    static_cast< Matrix<Scalar> >(*this) = other;
    this->_is_upper = other._is_upper;
    this->_is_unit_triangular = other._is_unit_triangular;

    return *this;
  }


  /**
   * Returns true if the matrix is a upper or lower triangular matrix.
   */
  inline bool isUpper() const
  {
    if (this->_transposed)
      return ! this->_is_upper;
    return this->_is_upper;
  }


  /**
   * Retruns true, if the matrix has a unit diagonal.
   */
  inline bool isUnit() const
  {
    return this->_is_unit_triangular;
  }


  /**
   * Constructs a weak view of this matrix in column-major storage order.
   */
  inline TriMatrix<Scalar> colMajor() const
  {
    // If storage order is column-major:
    if (! this->_is_rowmajor)
      return *this;

    return TriMatrix<Scalar>(Matrix<Scalar>::colMajor(),
                             ! this->_is_unit_triangular, this->_is_unit_triangular);
  }
};


/**
 * Constructs an upper-triangular matrix from the given general matrix.
 */
template <class Scalar>
inline TriMatrix<Scalar> triu(Matrix<Scalar> &matrix, bool unit=false)
{
  return TriMatrix<Scalar>(matrix, true, unit);
}


/**
 * Constructs a lower-triangular matrix from the given general matrix.
 */
template <class Scalar>
inline TriMatrix<Scalar> tril(Matrix<Scalar> &matrix, bool unit=false)
{
  return TriMatrix<Scalar>(matrix, false, unit);
}


}


#endif // __LINALG_TRIMATRIX_HH__
