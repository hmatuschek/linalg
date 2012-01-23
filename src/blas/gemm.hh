#ifndef __LINALG_BLAS_GEMM_HH__
#define __LINALG_BLAS_GEMM_HH__

extern "C" {
extern void dgemm_(const char *transa, const char *transb,
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
 * @ingroup linalg
 */
void gemm(const double alpha, const Matrix<double> &a, const Matrix<double> &b,
          const double beta, Matrix<double> &c)
{
  // Check if shape matches:
  if ( (a.cols() != b.rows()) || (a.rows() != c.rows()) || (b.cols() != c.cols()))
  {
    Linalg::IndexError err;
    err << "Shape mismatch.";
    throw err;
  }

  char transa='N';
  char transb='N';
  int M = a.rows();
  int N = b.cols();
  int K = a.cols();
  int lda = a.leading_dimension();
  int ldb = b.leading_dimension();
  int ldc = c.leading_dimension();

  // Perform operation:
  dgemm_(&transa, &transb, &M, &N, &K,
         &alpha, *a, &lda, *b, &ldb,
         &beta, *c, &ldc);

  // Done.
}


}
}

#endif
