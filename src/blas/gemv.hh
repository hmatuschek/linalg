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


#include "blas/utils.hh"
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
inline void gemv(const double alpha, const Matrix<double> &A, const Vector<double> &x,
          const double beta, Vector<double> &y)
{
  // Get matrix in column order (Fortran)
  char trans = 'N';
  Matrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, trans);

  LINALG_SHAPE_ASSERT(A.cols() == x.dim());
  LINALG_SHAPE_ASSERT(A.rows() == y.dim());

  int m      = Acol.rows();
  int n      = Acol.cols();
  int lda    = BLAS_LEADING_DIMENSION(Acol);
  int incx   = BLAS_INCREMENT(x);
  int incy   = BLAS_INCREMENT(y);

  dgemv_(&trans, &m, &n, &alpha, Acol.ptr(), &lda, x.ptr(), &incx, &beta, y.ptr(), &incy);
}


}
}
#endif // __LINALG_BLAS_GEMV_HH__
