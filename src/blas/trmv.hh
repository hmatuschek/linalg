/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_TRMV_HH__
#define __LINALG_BLAS_TRMV_HH__

#include "vector.hh"
#include "trimatrix.hh"
#include "blas/utils.hh"


// Iternface to Fortran function:
extern "C" {
void dtrmv_(char *uplo, char *transa, char *diag, int *n,
            const double *a, int *lda, double *x, int *incx);
}


namespace Linalg {
namespace Blas {


/**
 * Interfaces the DTRMV BLAS function.
 *
 * Calculates:
 * \f[x \leftarrow A\cdot x\f]
 * where @c A is a unit or non-unit, upper- or lower-triangular square matrix and @c x is a
 * vector.
 *
 * @param A Specifies the triangular matrix A.
 * @param upper If true, A is upper triangular, else A is lower triangular.
 * @param diag If true, A has unit diagonal.
 * @param x The vector x.
 * @throws ShapeError If A is not a square matrix, rows(A) != dim(x).
 *
 *
 * @ingroup blas2
 */
inline void trmv(const Matrix<double> &A, bool upper, bool diag, Vector<double> &x)
throw (ShapeError)
{
  // Assert Acol is square:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());
  LINALG_SHAPE_ASSERT(A.rows() == x.dim());

  Matrix<double> Acol = A;
  char transa = 'N';
  // Make A column-major:
  if (A.isRowMajor()) {
    upper = !upper;
    transa = 'T';
    Acol = Acol.t();
  }

  // Get flags and dimensions:
  char uplo  = upper ? 'U' : 'L';
  char d     = diag ? 'U' : 'N';
  int  N     = Acol.cols();
  int  lda   = Acol.strides()[1];
  int  incx  = BLAS_INCREMENT(x);

  // Call Fortran function
  dtrmv_(&uplo, &transa, &d, &N, Acol.ptr(), &lda, x.ptr(), &incx);
}


/**
 * Interfaces the DTRMV BLAS function.
 *
 * Calculates:
 * \f[x \leftarrow A\cdot x\f]
 * where @c A is a unit or non-unit, upper- or lower-triangular square matrix and @c x is a
 * vector.
 *
 * @ingroup blas2
 */
inline void trmv(const TriMatrix<double> &A, Vector<double> &x)
{
  trmv(A, A.isUpper(), A.hasUnitDiag(), x);
}


}
}
#endif // __LINALG_BLAS_TRMV_HH__
