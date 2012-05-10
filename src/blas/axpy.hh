#ifndef __LINALG_BLAS_AXPY_HH__
#define __LINALG_BLAS_AXPY_HH__

#include <vector.hh>
#include <simd.hh>


namespace Linalg {
namespace Blas {


/**
 * SIMD implementation of @c Linalg::Blas::axpy for dense vectors.
 *
 * @todo This function is untested.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline void __axpy_dense(const Scalar &alpha, size_t N, const Scalar *x, Scalar *y)
{
  size_t N_elm = SIMDTraits<Scalar>::num_elements;
  SIMDTraits<Scalar>::uvector alpha_vec;
  SIMDTraits<Scalar>::uvector *x_ptr = (SIMDTraits<Scalar>::uvector *)x;
  SIMDTraits<Scalar>::uvector *y_ptr = (SIMDTraits<Scalar>::uvector *)y;

  for (size_t i=0; i<N_elm; i++)
    alpha_vec.d[i] = alpha;

  size_t N_steps = N/N_elm;
  size_t N_rem   = N%N_elm;

  for (size_t i=0; i<N_steps; i++) {
    y_ptr->v += alpha_vec->v*x_ptr->v;
  }

  for (size_t i=0; i<N_rem; i++) {
    y_ptr->d[i] += alpha*x_ptr->d[i];
  }
}


/**
 * Direct implementation of @c Linalg::Blas::axpy for non-dense vectors.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline void __axpy_incremental(const Scalar &alpha, size_t N, const Scalar *x, size_t x_inc, Scalar *y, size_t y_inc)
{
  for (size_t i=0; i<N; i++, x+=x_inc, y += y_inc) {
    (*y) += alpha * (*x);
  }
}



/**
 * Computes \f$y' = \alpha x + y\f$ with x,y being vectors of same dimension.
 *
 * @ingroup blas1
 */
template<class Scalar>
inline void axpy(const Scalar &alpha, const Vector<Scalar> &x, Vector<Scalar> &y)
{
  LINALG_SHAPE_ASSERT(x.dim() == y.dim());

  if (0.0 == alpha)
    return;

  if ( (1 == x.strides()[0]) && (1 == y.strides()[0]))
    __axpy_dense(alpha, x.dim(), x.ptr(), y.ptr());
  else
    __axpy_incremental(alpha, x.dim(), x.ptr(), x.strides()[0], y.ptr(), y.strides()[0]);
}


}
}
#endif // __LINALG_BLAS_AXPY_HH__
