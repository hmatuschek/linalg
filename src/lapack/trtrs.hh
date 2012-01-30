/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_LAPACK_TRTRS_HH__
#define __LINALG_LAPACK_TRTRS_HH__

#include "blas/utils.hh"
#include "blas/trsm.hh"
#include "matrix.hh"
#include "trimatrix.hh"


namespace Linalg {
namespace Lapack {

/**
 * Solves the triangular system:
 *
 * \f[A\cdot X = B\f]
 *
 * Where @c A is a unit- or non-unit upper- or lower-triangular matrix and @c B is a general
 * matrix.
 *
 * The LAPACK ?TRTRS implementation was not flexible enougth, so I decided to reimplement
 * this function as a template. The LAPACK function simply checks if there is a 0 on the
 * diagonal of @c A. If not it calls BLASs ?TRSM, so do I.
 *
 * See also: @c Blas::trsm
 *
 * @throws Linalg::SingularMatrixError If matrix @c A is signular.
 * @throws Linalg::ShapeError If the shapes of @c A and @c B do not match.
 *
 * @ingroup lapack
 */
template <class Scalar>
inline void trtrs(const TriMatrix<Scalar> &A, Matrix<Scalar> &B)
throw (SingularMatrixError, ShapeError)
{
  // Ensure column-major representation of matrices:
  TriMatrix<Scalar> Acol = BLAS_TO_COLUMN_MAJOR(A);
  Matrix<Scalar> Bcol = BLAS_TO_COLUMN_MAJOR(B);

  // A must be quadratic:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());
  LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.rows());

  bool left = true;
  // If B is transposed, transpose A & B and swap sides:
  if (Bcol.isTransposed()) {
    Bcol = Bcol.t();
    Acol = Acol.t();
    left = !left;
  }

  // Check if A has zero on diagonal:
  if (! A.hasUnitDiag()) {
    for (size_t i=0; i<Acol.rows(); i++) {
      if (0 == Acol(i,i)) {
        SingularMatrixError err;
        err << "Signular matrix: " << i << "-th diagonal element of A is 0!";
        throw err;
      }
    }
  }

  // Just call trsm() ...
  Blas::trsm(Acol, 1.0, Bcol, left);
}


}
}
#endif // TRTRS_HH
