/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_TRIMATRIX_HH__
#define __LINALG_TRIMATRIX_HH__

#include "matrix.hh"

namespace Linalg {

/**
 * Implements a view on @c Matrix as a unit or non-unit upper- or lower-triangular matrix.
 *
 * @ingroup matrix
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
  TriMatrix(const Matrix<Scalar> &matrix, bool upper, bool unit)
    : Matrix<Scalar>(matrix), _is_upper(upper), _is_unit_triangular(unit)
  {
    // Pass...
  }


  /**
   * Copy constructor.
   */
  TriMatrix(const TriMatrix<Scalar> &other)
    : Matrix<Scalar>(other),
      _is_upper(other._is_upper), _is_unit_triangular(other._is_unit_triangular)
  {
    // Pass...
  }


  /**
   * Assignment operation.
   */
  TriMatrix<Scalar> &operator= (const TriMatrix<Scalar> &other)
  {
    // Assign Matrix:
    static_cast< Matrix<Scalar> &>(*this) = other;
    this->_is_upper = other._is_upper;
    this->_is_unit_triangular = other._is_unit_triangular;

    return *this;
  }


  Scalar operator() (size_t i, size_t j)
  {
    return Matrix<Scalar>::operator ()(i,j);
  }


  Scalar operator() (size_t i, size_t j) const
  {
    if (this->isUpper() && i <= j) {
      if (this->hasUnitDiag() && i==j) {
        return 1.0;
      }
      return Matrix<Scalar>::operator ()(i,j);
    } else if (!this->isUpper() && i >= j) {
      if (this->hasUnitDiag() && i==j) {
        return 1.0;
      }
      return Matrix<Scalar>::operator ()(i,j);
    }

    return 0.0;
  }


  /**
   * Returns true if the matrix is a upper or lower triangular matrix.
   */
  inline bool isUpper() const
  {
    return this->_is_upper;
  }


  /**
   * Retruns true, if the matrix has a unit diagonal.
   */
  inline bool hasUnitDiag() const
  {
    return this->_is_unit_triangular;
  }


  /**
   * Creates a @c TriMatrix view of this TriMatrix as a transposed.
   */
  inline TriMatrix<Scalar> t() const
  {
    return TriMatrix(Matrix<Scalar>::t(), !_is_upper, _is_unit_triangular);
  }
};


/**
 * Constructs an upper-triangular matrix from the given general matrix.
 */
template <class Scalar>
inline TriMatrix<Scalar> triu(const Matrix<Scalar> &matrix, bool unit=false)
{
  return TriMatrix<Scalar>(matrix, true, unit);
}


/**
 * Constructs a lower-triangular matrix from the given general matrix.
 */
template <class Scalar>
inline TriMatrix<Scalar> tril(const Matrix<Scalar> &matrix, bool unit=false)
{
  return TriMatrix<Scalar>(matrix, false, unit);
}


}


#endif // __LINALG_TRIMATRIX_HH__
