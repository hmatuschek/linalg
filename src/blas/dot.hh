/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_DOT_HH__
#define __LINALG_BLAS_DOT_HH__

#include "vector.hh"

extern "C" {
double ddot_(const int *N, const double *DX, const int *INCX, const double *DY, const int *INCY);
}


namespace Linalg {
namespace Blas {


/**
 * Calculates simple dot product:
 *
 * \f[dot(x,y) = x^T\cdot y\f]
 *
 * @ingroup blas1
 */
double dot(const Vector<double> &x, const Vector<double> &y) {
  if (x.dim() != y.dim()) {
    IndexError err; err << "Dimension mismatch."; throw err;
  }

  int N = x.dim();
  int INCX = x.stride();
  int INCY = y.stride();

  return ddot_(&N, *x, &INCX, *y, &INCY);
}


}
}
#endif // __LINALG_BLAS_DOT_HH__
