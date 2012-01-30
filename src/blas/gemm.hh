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
  // Get matrices in column-major from:
  Matrix<double> Acol = BLAS_TO_COLUMN_MAJOR(A);
  Matrix<double> Bcol = BLAS_TO_COLUMN_MAJOR(B);
  Matrix<double> Ccol = BLAS_TO_COLUMN_MAJOR(C);

  // If C is transposed:
  if (Ccol.isTransposed()){
    // Swap A & B
    Matrix<double> tmp(Acol.t());
    Acol = Bcol.t(); Bcol = tmp;
    // Transpose C:
    Ccol = Ccol.t();
  }

  LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.rows());
  LINALG_SHAPE_ASSERT(Acol.rows() == Ccol.rows());
  LINALG_SHAPE_ASSERT(Bcol.cols() == Ccol.cols());

  char transa = BLAS_TRANSPOSED_FLAG(Acol);
  char transb = BLAS_TRANSPOSED_FLAG(Bcol);
  int M       = Acol.rows();
  int N       = Bcol.cols();
  int K       = Acol.cols();
  int lda     = BLAS_LEADING_DIMENSION(Acol);
  int ldb     = BLAS_LEADING_DIMENSION(Bcol);
  int ldc     = BLAS_LEADING_DIMENSION(Ccol);

  // Perform operation:
  dgemm_(&transa, &transb, &M, &N, &K,
         &alpha, *Acol, &lda,
         *Bcol, &ldb,
         &beta, *Ccol, &ldc);

  // Done.
}


}
}

#endif
