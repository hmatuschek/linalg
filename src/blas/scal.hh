/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_SCAL_HH__
#define __LINALG_BLAS_SCAL_HH__

#include "vector.hh"
#include "utils.hh"
#include "simd.hh"


namespace Linalg {
namespace Blas {


/**
 * Optimized internal function to scale a vector using SIMD instructions if the vector is dense
 * (increment 1).
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline void __scal_dense(size_t N, const Scalar &a, Scalar *x) {
  typename SIMDTraits<Scalar>::uvector a_vec;
  typename SIMDTraits<Scalar>::uvector *x_ptr = (typename SIMDTraits<Scalar>::uvector *)x;

  size_t N_elm  = SIMDTraits<Scalar>::num_elements;
  size_t N_step = N/N_elm;
  size_t N_rem  = N%N_elm;

  // Initialize constant vector a_vec
  for (size_t i=0; i<N_elm; i++) {
    a_vec.d[i] = a;
  }

  // Perform operations on vectors:
  for (size_t i=0; i<N_step; i++, x_ptr++) {
    x_ptr->v *= a_vec.v;
  }

  // Handle remaining elements:
  for (size_t i=0; i<N_rem; i++) {
    x_ptr->d[i] *= a;
  }
}


/**
 * Internal function to scale a vector.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline void __scal_incremental(size_t N, const Scalar &a, Scalar *x, size_t xinc) {
  for (size_t i=0; i<N; i++, x += xinc) {
    (*x) *= a;
  }
}


/**
 * Scales a vector inplace.
 *
 * \f[x' = ax\f]
 *
 * This function uses @c __scal_dense if x is dense (increment 1) and @c __scal_incremental
 * otherwise.
 *
 * @ingroup blas1
 */
template <class Scalar>
inline void scal(const Scalar &a, Vector<Scalar> &x) {
  size_t N    = x.dim();
  size_t xinc = BLAS_INCREMENT(x);

  if (1 == xinc)
    __scal_dense(N, a, x.ptr());
  else
    __scal_incremental(N, a, x.ptr(), xinc);
}


}
}
#endif // ____LINALG_BLAS_SCAL_HH__
