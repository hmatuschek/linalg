#ifndef __LINALG_BLAS_GEMM_HH__
#define __LINALG_BLAS_GEMM_HH__

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
 * @bug This function did not pass the unit-tests.
 *
 * @ingroup blas3
 */
void gemm(double alpha, const Matrix<double> &A, const Matrix<double> &B,
          double beta, Matrix<double> &C)
{
  // Get matrices in column-major from:
  Matrix<double> Acol = A.colMajor();
  Matrix<double> Bcol = B.colMajor();
  Matrix<double> Ccol = C.colMajor();

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

  char transa='N';
  char transb='N';

  if (Acol.isTransposed())
    transa = 'T';
  if (Bcol.isTransposed())
    transb = 'T';

  int M = Acol.rows();
  int N = Bcol.cols();
  int K = Acol.cols();
  int lda = Acol.stride();
  int ldb = Bcol.stride();
  int ldc = Ccol.stride();

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
