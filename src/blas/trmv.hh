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
 * @ingroup blas2
 */
inline void trmv(const TriMatrix<double> &A, Vector<double> &x)
{
  // Assert Acol is square:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());
  LINALG_SHAPE_ASSERT(A.rows() == x.dim());

  // Make A column-major:
  char transa = 'N';
  TriMatrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, transa);

  // Get flags and dimensions:
  char uplo  = BLAS_UPLO_FLAG(Acol);
  char diag  = BLAS_UNIT_DIAG_FLAG(Acol);
  int  N     = BLAS_NUM_COLS(Acol, transa);
  int  lda   = BLAS_LEADING_DIMENSION(Acol);
  int  incx  = BLAS_INCREMENT(x);

  // Call Fortran function
  dtrmv_(&uplo, &transa, &diag, &N, Acol.ptr(), &lda, x.ptr(), &incx);
}


}
}
#endif // __LINALG_BLAS_TRMV_HH__
