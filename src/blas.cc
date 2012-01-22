#include "blas.hh"
#include "exception.hh"

#include <iostream>

/*
 * C interface to FORTRAN functions (direct access).
 */
extern "C" {
extern void dgemv_(const char *trans, const int *m, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx,
                   const double *beta, double *y, const int *incy);

extern void dgemm_(const char *transa, const char *transb,
                   const int *m, const int *n, const int *k,
                   const double *alpha, const double *a, const int *lda,
                   const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc);
}


using namespace Fluc;
using namespace Fluc::Linalg;



void
Linalg::gemv(const double alpha, const Matrix<double> &a, const Vector<double> &x,
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
  int incx = 1;
  int incy = 1;

  dgemv_(&trans, &m, &n, &alpha, *a, &lda, *x, &incx, &beta, *y, &incy);
}


void
Linalg::gemm(const double alpha, const Matrix<double> &a, const Matrix<double> &b,
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



Matrix<double>
Linalg::dot(const Matrix<double> &a, const Matrix<double> &b)
{
  // Allocate matrix for C = A*B:
  Matrix<double> c(a.rows(), b.cols());

  // Perform in-place:
  gemm(1., a, b, 0., c);

  // Done.
  return c;
}


Vector<double>
Linalg::dot(const Matrix<double> &a, const Vector<double> &b)
{
  // Allocate
  Vector<double> c(a.rows());

  // perform operation in-place:
  gemv(1, a, b, 0., c);

  // Done.
  return c;
}


