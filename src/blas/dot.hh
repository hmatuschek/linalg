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
 * @ingroup blas_internal
 */
template <class Scalar>
Scalar __dot_dense(size_t N, const Scalar *x, const Scalar *y)
{
  typename SIMDTraits<Scalar>::uvector res;
  const typename SIMDTraits<Scalar>::uvector *x_ptr = (const typename SIMDTraits<Scalar>::uvector *)x;
  const typename SIMDTraits<Scalar>::uvector *y_ptr = (const typename SIMDTraits<Scalar>::uvector *)y;

  // get n-blocks, and remainder
  size_t N_elm  = SIMDTraits<Scalar>::num_elements;
  size_t N_step = N/N_elm;
  size_t N_rem  = N%N_elm;

  // Initialize result vector:
  for (size_t i=0; i<N_elm; i++) {
    res.d[i] = Scalar(0);
  }

  // Work on SIMD vectors
  for (size_t i=0; i<N_step; i++, x_ptr++, y_ptr++) {
    res.v += x_ptr->v * y_ptr->v;
  }

  // Handle last elements...
  for (size_t i=0; i<N_rem; i++) {
    res.d[i] += x_ptr->d[i] * y_ptr->d[i];
  }

  // Compute sum:
  for (size_t i=1; i<N_elm; i++) {
    res.d[0] += res.d[i];
  }

  return res.d[0];
}


/**
 * Internal used function, calculateing the inner product as \f$x^Ty\f$.
 *
 * @note This function does no dimension check on x & y.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline Scalar __dot_incremental(size_t N, const Scalar *x, size_t inc_x, const Scalar *y, size_t inc_y)
{
  Scalar r = Scalar(0);

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
    return __dot_dense<Scalar>(N, x.ptr(), y.ptr());
  }

  // otherwise use incremental operation
  return __dot_incremental<Scalar>(N, x.ptr(), incx, y.ptr(), incy);
}


}
}

#endif // __LINALG_BLAS_DOT_HH__
