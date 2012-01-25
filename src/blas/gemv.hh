/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_GEMV_HH__
#define __LINALG_BLAS_GEMV_HH__

extern "C" {
extern void dgemv_(const char *trans, const int *m, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx,
                   const double *beta, double *y, const int *incy);
}

#include "matrix.hh"


namespace Linalg {
namespace Blas {


/**
 * Direct access to the dgemv BLAS function.
 *
 * Calculates in-place:
 * \f[y = \alpha*A*x + \beta*y\f]
 *
 * @ingroup blas2
 */
void gemv(const double alpha, const Matrix<double> &A, const Vector<double> &x,
          const double beta, Vector<double> &y)
{
  // Get matrix in column order (Fortran)
  Matrix<double> Acol(A.colMajor());

  LINALG_SHAPE_ASSERT(Acol.cols() == x.dim());
  LINALG_SHAPE_ASSERT(Acol.rows() == y.dim());

  char trans = 'N';
  int m = Acol.rows();
  int n = Acol.cols();

  if (Acol.isTransposed()) {
    trans = 'T';
    std::swap(m,n);
  }

  int lda = Acol.stride();
  int incx = x.stride();
  int incy = y.stride();

  dgemv_(&trans, &m, &n, &alpha, *Acol, &lda, *x, &incx, &beta, *y, &incy);
}


}
}
#endif // __LINALG_BLAS_GEMV_HH__
