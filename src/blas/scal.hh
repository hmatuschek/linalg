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
 * @ingroup blas1_internal
 */
inline void __scal_dense(size_t N, double a, double *x) {
  __vec2df *x_ptr = (__vec2df *)x;
  __v2df a_vec; a_vec = (__v2df) {a, a};

  size_t N2 = N/2;
  size_t r  = N%2;

  for (size_t i=0; i<N2; i++, x_ptr++) {
    x_ptr->v *= a_vec;
  }

  if (1 == r) {
    x_ptr->d[0] *= a;
  }
}


/**
 * Internal function to scale a vector.
 *
 * @ingroup blas1_internal
 */
inline void __scal_incremental(size_t N, const double &a, double *x, size_t xinc) {
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
inline void scal(const double &a, Vector<double> &x) {
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
