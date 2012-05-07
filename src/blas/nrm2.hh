/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_NRM2_HH__
#define __LINALG_BLAS_NRM2_HH__


#include "blas/utils.hh"
#include "vector.hh"
#include "simd.hh"
#include <cmath>


namespace Linalg {
namespace Blas {


/**
 * Internal function to perform @c nrm2 using SIMD instructions.
 *
 * @ingroup blas_internal
 */
inline double __nrm2_dense(size_t N, const double *x) {
  __vec2df res; res.v = (__v2df){0.0, 0.0};
  const __vec2df *x_ptr = (const __vec2df *)x;

  size_t N2 = N/2;
  size_t r  = N%2;

  for(size_t i=0; i<N2; i++, x_ptr++) {
    res.v += x_ptr->v * x_ptr->v;
  }

  if (r) {
    res.d[0] += x_ptr->d[0] * x_ptr->d[0];
  }

  res.d[0] += res.d[1];
  return sqrt(res.d[0]);
}


/**
 * Internal function to perform @c nrm2 for non-dense vectors.
 *
 * @ingroup blas_internal
 */
inline double __nrm2_incremental(size_t N, const double *x, size_t incx) {
  double res = 0.0;

  for (size_t i=0; i<N; i++, x+=incx) {
    res += (*x) * (*x);
  }

  return sqrt(res);
}



/**
 * Calculates 2-norm of a vector.
 *
 * \f[nrm2(x) = \sqrt{dot(x,x)} \f]
 *
 * @ingroup blas1
 */
inline double nrm2(const Vector<double> &x)
{
  int N   = BLAS_DIMENSION(x);
  int INC = BLAS_INCREMENT(x);

  if (1 == INC) {
    return __nrm2_dense(N, x.ptr());
  }

  return __nrm2_incremental(N, x.ptr(), INC);
}


}
}
#endif // __LINALG_BLAS_NRM2_HH__
