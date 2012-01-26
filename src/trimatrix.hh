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
    if (matrix.isTransposed())
      this->_is_upper = ! this->_is_upper;
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
  inline bool hasUnitDiag() const
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

    TriMatrix<Scalar> col(Matrix<Scalar>::colMajor(),
                          this->isUpper(), this->hasUnitDiag());
    col._is_upper = ! this->_is_upper;
    return col;
  }


  /**
   * Creates a @c TriMatrix view of this TriMatrix as a transposed.
   */
  inline TriMatrix<Scalar> t() const
  {
    TriMatrix<Scalar> trans(*this); trans._transposed = !trans._transposed;
    return trans;
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
