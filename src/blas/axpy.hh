#ifndef __LINALG_BLAS_AXPY_HH__
#define __LINALG_BLAS_AXPY_HH__

#include <vector.hh>
#include <simd.hh>


namespace Linalg {
namespace Blas {


/**
 * SIMD implementation of @c Linalg::Blas::axpy for dense vectors.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline void __axpy_dense(const Scalar &alpha, size_t N, const Scalar *x, Scalar *y)
throw ()
{
  typename SIMDTraits<Scalar>::uvector alpha_vec;
  const typename SIMDTraits<Scalar>::uvector *x_ptr = (const typename SIMDTraits<Scalar>::uvector *)x;
  typename SIMDTraits<Scalar>::uvector *y_ptr = (typename SIMDTraits<Scalar>::uvector *)y;

  size_t N_elm = SIMDTraits<Scalar>::num_elements;
  size_t N_steps = N/N_elm;
  size_t N_rem   = N%N_elm;

  // Initialize alpha_vector
  for (size_t i=0; i<N_elm; i++)
    alpha_vec.d[i] = alpha;

  // Perform on vectors:
  for (size_t i=0; i<N_steps; i++, x_ptr++, y_ptr++)
    y_ptr->v += alpha_vec.v * x_ptr->v;

  // Perform on remaining elements
  for (size_t i=0; i<N_rem; i++)
    y_ptr->d[i] += alpha*x_ptr->d[i];
}


/**
 * Direct implementation of @c Linalg::Blas::axpy for non-dense vectors.
 *
 * @ingroup blas_internal
 */
template <class Scalar>
inline void
__axpy_incremental(const Scalar &alpha, size_t N, const Scalar *x, size_t x_inc,
                   Scalar *y, size_t y_inc) throw ()
{
  for (size_t i=0; i<N; i++, x+=x_inc, y+=y_inc) {
    (*y) += alpha * (*x);
  }
}


/**
 * Computes \f$y' = \alpha x + y\f$ with x,y being vectors of same dimension.
 *
 * This function calls @c Linalg::Blas::__axpy_dense if both vectors are stored densly, otherwise
 * @c Linalg::Blas::__axpy_incremental is called.
 *
 * @note In contrast to @c Linalg::Blas::axpy, this function does no dimension check on the vectors.
 *
 * @ingroup blas_internal
 */
template<class Scalar>
inline void __axpy(const Scalar &alpha, const Vector<Scalar> &x, Vector<Scalar> &y)
throw ()
{
  if (0.0 == alpha)
    return;

  if ( (1 == x.strides()[0]) && (1 == y.strides()[0]))
    __axpy_dense(alpha, x.dim(), x.ptr(), y.ptr());
  else
    __axpy_incremental(alpha, x.dim(), x.ptr(), x.strides()[0], y.ptr(), y.strides()[0]);
}



/**
 * Computes \f$y' = \alpha x + y\f$ with x,y being vectors of same dimension.
 *
 * @throws ShapeError If dim(x) != dim(y).
 *
 * @ingroup blas1
 */
template<class Scalar>
inline void axpy(const Scalar &alpha, const Vector<Scalar> &x, Vector<Scalar> &y)
throw (ShapeError)
{
  LINALG_SHAPE_ASSERT(x.dim() == y.dim());

  __axpy(alpha, x, y);
}


}
}
#endif // __LINALG_BLAS_AXPY_HH__
