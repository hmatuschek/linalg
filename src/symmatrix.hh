/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_SYMMATRIX_HH__
#define __LINALG_SYMMATRIX_HH__

#include "matrix.hh"

namespace Linalg {

/**
 * Implements a view on a @c Matrix as a symtric matrix, where only the upper- or lower-triangular
 * part is used.
 *
 * The structure of this matrix is identical to a @c TriMatrix, but it has some different
 * semantic, therefore there are two classes.
 *
 * @ingroup matrix
 */
template <class Scalar>
class SymMatrix : public Matrix<Scalar>
{
protected:
  /**
   * Specifies if 1/2 of the matrix is storred in the upper or lower half.
   */
  bool _is_upper;


public:
  /**
   * Constructs a triangular matrix-view from the given marix.
   */
  SymMatrix(const Matrix<Scalar> &matrix, bool upper)
    : Matrix<Scalar>(matrix), _is_upper(upper)
  {
    if (matrix.isTransposed())
      this->_is_upper = ! this->_is_upper;
  }


  /**
   * Copy constructor.
   */
  SymMatrix(const SymMatrix<Scalar> &other)
    : Matrix<Scalar>(other), _is_upper(other._is_upper)
  {
    // Pass...
  }


  /**
   * Assignment operation.
   */
  SymMatrix<Scalar> &operator= (const SymMatrix<Scalar> &other)
  {
    // Assign Matrix:
    static_cast< Matrix<Scalar> &>(*this) = other;
    this->_is_upper = other._is_upper;

    return *this;
  }


  /**
   * Returns true if 1/2 of the matrix is stored as a upper or lower triangular matrix.
   */
  inline bool isUpper() const
  {
    if (this->_transposed)
      return ! this->_is_upper;
    return this->_is_upper;
  }


  /**
   * Constructs a weak view of this matrix in column-major storage order.
   */
  inline SymMatrix<Scalar> colMajor() const
  {
    // If storage order is column-major:
    if (! this->_is_rowmajor)
      return *this;

    SymMatrix<Scalar> col(Matrix<Scalar>::colMajor(),
                          this->isUpper());
    col._is_upper = ! this->_is_upper;
    return col;
  }


  /**
   * Creates a @c SymMatrix view of this SymMatrix as a transposed.
   */
  inline TriMatrix<Scalar> t() const
  {
    SymMatrix<Scalar> trans(*this); trans._transposed = !trans._transposed;
    return trans;
  }
};


/**
 * Constructs a symetric matrix using the upper-triangular part from the given general matrix.
 */
template <class Scalar>
inline SymMatrix<Scalar> symu(const Matrix<Scalar> &matrix)
{
  return SymMatrix<Scalar>(matrix, true);
}


/**
 * Constructs a symetric matrix using the lower-triangular part from the given general matrix.
 */
template <class Scalar>
inline SymMatrix<Scalar> syml(const Matrix<Scalar> &matrix)
{
  return SymMatrix<Scalar>(matrix, false);
}


}


#endif // __LINALG_TRIMATRIX_HH__
