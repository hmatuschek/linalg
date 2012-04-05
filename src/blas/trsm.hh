/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_TRSM_HH__
#define __LINALG_BLAS_TRSM_HH__


extern "C" {
void dtrsm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG,
            const int *M, const int *N,
            const double *alpha, const double *A, const int *LDA,
            double *B, const int *LDB);
}


#include "blas/utils.hh"
#include "matrix.hh"
#include "trimatrix.hh"


namespace Linalg {
namespace Blas{

/**
 * Interface to BLAS level 3 DTRSM function.
 *
 * Solves the triangular system \f$ A\cdot X = \alpha B\f$ if @c left=true
 * or \f$X \cdot A = \alpha B\f$ if @c left=false.
 *
 * @ingroup blas3
 */
inline void
trsm(const TriMatrix<double> &A, double alpha, Matrix<double> &B, bool left=true)
throw (ShapeError)
{
  // Check shapes:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());
  if (left) {
    LINALG_SHAPE_ASSERT(A.cols() == B.rows())
  } else {
    LINALG_SHAPE_ASSERT(B.cols() == A.rows());
  }

  // Ensure A & B are column-major:
  char transa='N', transb='N';
  TriMatrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, transa);
  Matrix<double>    Bcol = B; BLAS_ENSURE_COLUMN_MAJOR(Bcol, transb);

  // If B is transposed -> transpose A & B and swap side:
  char uplo = BLAS_UPLO_FLAG(Acol);
  if ('T'==transb) {
    transa = BLAS_TRANSPOSE(transa);
    transb = BLAS_TRANSPOSE(transb);
    uplo   = BLAS_TRANSPOSE_UPLO(transa, uplo);
    left = !left;
  }

  // Get flags:
  char side   = left ? 'L' : 'R';
  char diag   = BLAS_UNIT_DIAG_FLAG(Acol);
  int  M      = BLAS_NUM_ROWS(Bcol, transb);
  int  N      = BLAS_NUM_COLS(Bcol, transb);
  int  lda    = BLAS_LEADING_DIMENSION(Acol);
  int  ldb    = BLAS_LEADING_DIMENSION(Bcol);

  // Done...
  dtrsm_(&side, &uplo, &transa, &diag, &M, &N, &alpha, Acol.ptr(), &lda, Bcol.ptr(), &ldb);
}


}
}

#endif // TRSM_HH
