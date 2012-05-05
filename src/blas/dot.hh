/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_DOT_HH__
#define __LINALG_BLAS_DOT_HH__

#include "blas/utils.hh"
#include "vector.hh"
#include "simd.hh"


namespace Linalg {
namespace Blas {


/**
 * Specialized, internal function to calculate the inner product of two dense vectors in an
 * efficient way using SIMD instructions.
 *
 * @note This function does no dimension checks on x and y.
 *
 * @ingroup blas1_internal
 */
inline double __dot_dense(size_t N, const double *x, const double *y)
{
  // Allocate and initialize result vector:
  __vec2df res; res.v = (__v2df){0.0, 0.0};
  // Cast pointer to x and y
  const __vec2df *x_ptr = (const __vec2df *)x;
  const __vec2df *y_ptr = (const __vec2df *)y;

  // get n-blocks, and remainder
  size_t N2 = N/2;
  size_t r  = N%2;

  // Work on two-element vectors
  for (size_t i=0; i<N2; i++, x_ptr++, y_ptr++) {
    res.v += x_ptr->v * y_ptr->v;
  }

  // Handle last element, if N is odd
  if (r) {
    res.d[0] += x_ptr->d[0] * y_ptr->d[0];
  }

  // Compute sum:
  res.d[0] += res.d[1];
  return res.d[0];
}


/**
 * Specialized, internal function to calculate the inner product of 2 dense vectors in an
 * efficient way using SIMD instructions.
 *
 * @note This function does no dimension checks on x and y.
 *
 * @ingroup blas1_internal
 */
inline float __dot_dense(size_t N, const float *x, const float *y)
{
  // Allocate and initialize result vector:
  __vec4sf res; res.v = (__v4sf){0.0f, 0.0f, 0.0f, 0.0f};
  // Cast pointer to x and y
  const __vec4sf *x_ptr = (const __vec4sf *)x;
  const __vec4sf *y_ptr = (const __vec4sf *)y;

  // get n-blocks, and remainder
  size_t N4 = N/4;
  size_t r  = N%4;

  // Work on two-element vectors
  for (size_t i=0; i<N4; i++, x_ptr++, y_ptr++) {
    res.v += x_ptr->v * y_ptr->v;
  }

  // Handle last elements
  for (size_t i=0; i<r; i++) {
    res.d[0] += x_ptr->d[i] * y_ptr->d[i];
  }

  // Compute sum:
  res.d[0] += res.d[1] + res.d[2] + res.d[3];
  return res.d[0];
}


/**
 * Internal used function, calculateing the inner product as \f$x^Ty\f$.
 *
 * @note This function does no dimension check on x & y.
 *
 * @ingroup blas1_internal
 */
inline double __dot_incremental(size_t N, const double *x, size_t inc_x, const double *y, size_t inc_y)
{
  double r = 0.0;

  for (size_t i=0; i<N; i++, x+=inc_x, y+=inc_y) {
    r +=  (*x) * (*y);
  }
  return r;
}


/**
 * Internal used function, calculateing the inner product as \f$x^Ty\f$.
 *
 * @note This function does no dimension check on x & y.
 *
 * @ingroup blas1_internal
 */
inline float __dot_incremental(size_t N, const float *x, size_t inc_x, const float *y, size_t inc_y)
{
  float r = 0;

  for (size_t i=0; i<N; i++, x+=inc_x, y+=inc_y) {
    r +=  (*x) * (*y);
  }
  return r;
}



/**
 * Calculates the dot product of two vectors x and y.
 *
 * \f[dot(x,y) = x^T\cdot y\f]
 *
 * If the vectors x and y are dense (increment 1), the optimized method @c __dot_dense is used
 * otherwise the method @c __dot_incremental is used.
 *
 * @ingroup blas1
 */
template <class Scalar>
inline Scalar dot(const Vector<Scalar> &x, const Vector<Scalar> &y)
{
  LINALG_SHAPE_ASSERT(x.dim() == y.dim());

  int N    = BLAS_DIMENSION(x);
  int incx = BLAS_INCREMENT(x);
  int incy = BLAS_INCREMENT(y);

  // If x and y are dense, use optimized methods:
  if ( (1 == incx) && (1 == incy) ) {
    return __dot_dense(N, x.ptr(), y.ptr());
  }

  // otherwise use incremental operation
  return __dot_incremental(N, x.ptr(), incx, y.ptr(), incy);
}


}
}

#endif // __LINALG_BLAS_DOT_HH__
