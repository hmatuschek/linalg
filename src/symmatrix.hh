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
#include "trimatrix.hh"


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
class SymMatrix : public TriMatrix<Scalar>
{
public:
  /**
   * Constructs a triangular matrix-view from the given marix.
   */
  SymMatrix(const Matrix<Scalar> &matrix, bool upper) : TriMatrix<Scalar>(matrix, upper, false) { }

  SymMatrix(const TriMatrix<Scalar> &other) : TriMatrix<Scalar>(other) { }

  /**
   * Copy constructor.
   */
  SymMatrix(const SymMatrix<Scalar> &other) : TriMatrix<Scalar>(other) { }

  inline SymMatrix<Scalar> &operator= (const SymMatrix &other)
  {
    TriMatrix<Scalar>::operator =(other);
    return *this;
  }

  inline SymMatrix<Scalar> t() const {
    return TriMatrix<Scalar>::t();
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
