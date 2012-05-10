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
 * Internal function to perform @c nrm2sq using SIMD instructions.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline Scalar __nrm2sq_dense(size_t N, const Scalar *x) {
  typename SIMDTraits<Scalar>::uvector res;
  const typename SIMDTraits<Scalar>::uvector *x_ptr = (const typename SIMDTraits<Scalar>::uvector *)x;

  size_t N_elm  = SIMDTraits<Scalar>::num_elements;
  size_t N_step = N/N_elm;
  size_t N_rem  = N%N_elm;

  // Initialize result vector:
  for (size_t i=0; i<N_elm; i++) {
    res.d[i] = Scalar(0);
  }

  // Perform operations on vectors:
  for(size_t i=0; i<N_step; i++, x_ptr++) {
    res.v += x_ptr->v * x_ptr->v;
  }

  // Handle remaining elements:
  for (size_t i=0; i<N_rem; i++) {
    res.d[i] += x_ptr->d[i] * x_ptr->d[i];
  }

  // calc result
  for (size_t i=1; i<N_elm; i++) {
    res.d[0] += res.d[i];
  }

  return res.d[0];
}


/**
 * Internal function to perform @c nrm2sq for non-dense vectors.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline Scalar __nrm2sq_incremental(size_t N, const Scalar *x, size_t incx) {
  Scalar res = Scalar(0);

  for (size_t i=0; i<N; i++, x+=incx) {
    res += (*x) * (*x);
  }

  return res;
}



/**
 * Calculates squared 2-norm of a vector.
 *
 * \f[nrm2sq(x) = dot(x,x) \f]
 *
 * @ingroup blas1
 */
template <class Scalar>
inline Scalar nrm2sq(const Vector<Scalar> &x)
{
  int N   = BLAS_DIMENSION(x);
  int INC = BLAS_INCREMENT(x);

  if (1 == INC) {
    return __nrm2sq_dense(N, x.ptr());
  }

  return __nrm2sq_incremental(N, x.ptr(), INC);
}



/**
 * Calculates 2-norm of a vector.
 *
 * \f[nrm2(x) = \sqrt{dot(x,x)} \f]
 *
 * @ingroup blas1
 */
template <class Scalar>
inline Scalar nrm2(const Vector<Scalar> &x)
{
  int N   = BLAS_DIMENSION(x);
  int INC = BLAS_INCREMENT(x);

  if (1 == INC) {
    return sqrt(__nrm2sq_dense(N, x.ptr()));
  }

  return sqrt(__nrm2sq_incremental(N, x.ptr(), INC));
}


}
}
#endif // __LINALG_BLAS_NRM2_HH__
