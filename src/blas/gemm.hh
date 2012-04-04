/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_GEMM_HH__
#define __LINALG_BLAS_GEMM_HH__


#include "blas/utils.hh"
#include "matrix.hh"


extern "C" {
void dgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);
}



namespace Linalg {
namespace Blas {


/**
 * Direct access to the dgemm LAPACK function.
 *
 * Calculates:
 * \f[ C = \alpha A * B + \beta * C \f]
 *
 * @ingroup blas3
 */
inline void gemm(double alpha, const Matrix<double> &A, const Matrix<double> &B,
          double beta, Matrix<double> &C)
{
  // Check matrix shape:
  LINALG_SHAPE_ASSERT(A.cols() == B.rows());
  LINALG_SHAPE_ASSERT(A.rows() == C.rows());
  LINALG_SHAPE_ASSERT(B.cols() == C.cols());

  // Get matrices in column-major from:
  char transa='N', transb='N', transc = 'N';
  Matrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, transa);
  Matrix<double> Bcol = B; BLAS_ENSURE_COLUMN_MAJOR(Bcol, transb);
  Matrix<double> Ccol = C; BLAS_ENSURE_COLUMN_MAJOR(Ccol, transc);

  if ('T' == transc) {
    std::swap(Acol, Bcol);
    std::swap<>(transa, transb);
    transa = BLAS_TRANSPOSE(transa);
    transb = BLAS_TRANSPOSE(transb);
    transc = BLAS_TRANSPOSE(transc);
  }

  int M       = BLAS_NUM_ROWS(Acol, transa);
  int N       = BLAS_NUM_COLS(Bcol, transb);
  int K       = BLAS_NUM_COLS(Acol, transa);
  int lda     = BLAS_LEADING_DIMENSION(Acol);
  int ldb     = BLAS_LEADING_DIMENSION(Bcol);
  int ldc     = BLAS_LEADING_DIMENSION(Ccol);

  // Perform operation:
  dgemm_(&transa, &transb, &M, &N, &K,
         &alpha, Acol.ptr(), &lda,
         Bcol.ptr(), &ldb,
         &beta, Ccol.ptr(), &ldc);

  // Done.
}


}
}

#endif
