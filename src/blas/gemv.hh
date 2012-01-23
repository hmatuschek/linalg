#ifndef __LINALG_BLAS_GEMV_HH__
#define __LINALG_BLAS_GEMV_HH__

extern "C" {
extern void dgemv_(const char *trans, const int *m, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx,
                   const double *beta, double *y, const int *incy);
}


namespace Linalg {
namespace Blas {


/**
 * Direct access to the dgemv LAPACK function.
 *
 * @ingroup linalg
 */
void gemv(const double alpha, const Matrix<double> &a, const Vector<double> &x,
          const double beta, Vector<double> &y)
{
  // Check if shape matches:
  if ((a.cols() != x.dim()) || (a.rows() != y.dim()))
  {
    Linalg::IndexError err;
    err << "Shape mismatch!";
    throw err;
  }

  char trans = 'N';
  int m = a.rows();
  int n = a.cols();
  int lda = a.leading_dimension();

  int incx = x.stride();
  int incy = y.stride();

  dgemv_(&trans, &m, &n, &alpha, *a, &lda, *x, &incx, &beta, *y, &incy);
}


}
}
#endif // __LINALG_BLAS_GEMV_HH__
