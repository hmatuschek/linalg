/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_TRMM_HH__
#define __LINALG_BLAS_TRMM_HH__

// Iternface to Fortran function:
extern "C" {
void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double *a, int *lda, double *b, int *ldb);
}



#include "trimatrix.hh"
#include "blas/utils.hh"


namespace Linalg {
namespace Blas {

/**
 * Implements an interface to the DTRMM function of BLAS.
 *
 * Calculates
 * \f[B = A\cdot B\f]
 * if @c left is true and
 * \f[B = B\cdot A\f]
 * if @c left is false.
 *
 * @todo This function is not tested yet.
 *
 * @param left Specifies if @c A is multiplied from left or right to @c B.
 * @param A Specifies a unit or non-unit, upper- or lower-triangular matrix.
 * @param B Specifies some general matrix.
 *
 * @ingroup blas3
 */
void trmm(bool left, double alpha, const TriMatrix<double> &A, Matrix<double> &B)
{
  // Make sure, A & B are in column-major form:
  TriMatrix<double> Acol(BLAS_TO_COLUMN_MAJOR(A));
  Matrix<double>    Bcol(BLAS_TO_COLUMN_MAJOR(B));

  // If Bcol is transposed: swap side and transpose Acol and Bcol:
  if (BLAS_IS_TRANSPOSED(Bcol)) {
    Acol = Acol.t();
    Bcol = Bcol.t();
    left = !left;
  }

  // Check dimensions:
  if (left) {
    LINALG_SHAPE_ASSERT(Acol.rows() == Bcol.rows());
    LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.rows());
  } else {
    LINALG_SHAPE_ASSERT(Acol.rows() == Bcol.cols());
    LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.cols());
  }

  char SIDE    = left ? 'L' : 'R';
  char UPLO    = BLAS_UPLO_FLAG(Acol);
  char TRANSA  = BLAS_TRANSPOSED_FLAG(Acol);
  char DIAG    = BLAS_UNIT_DIAG_FLAG(Acol);
  int  M       = BLAS_NUM_ROWS(Bcol);
  int  N       = BLAS_NUM_COLS(Bcol);
  double ALPHA = alpha;
  int LDA      = BLAS_LEADING_DIMENSION(Acol);
  int LDB      = BLAS_LEADING_DIMENSION(Bcol);

  // Call fortran function:
  dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, *Acol, &LDA, *Bcol, &LDB);
}

}
}

#endif // __LINALG_BLAS_TRMM_HH__
