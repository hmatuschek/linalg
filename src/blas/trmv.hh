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
void trmv(const TriMatrix<double> &A, Vector<double> &x)
{
  // Make A column-major:
  const TriMatrix<double> Acol = BLAS_TO_COLUMN_MAJOR(A);

  // Assert Acol is square:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());
  LINALG_SHAPE_ASSERT(Acol.rows() == x.dim());

  // Get flags and dimensions:
  char UPLO  = BLAS_UPLO_FLAG(Acol);
  char TRANS = BLAS_TRANSPOSED_FLAG(Acol);
  char DIAG  = BLAS_UNIT_DIAG_FLAG(Acol);
  int  N     = BLAS_NUM_COLS(Acol);
  int  LDA   = BLAS_LEADING_DIMENSION(Acol);
  int  INCX  = BLAS_INCREMENT(x);

  // Call Fortran function
  dtrmv_(&UPLO, &TRANS, &DIAG, &N, *Acol, &LDA, *x, &INCX);
}


}
}
#endif // __LINALG_BLAS_TRMV_HH__
