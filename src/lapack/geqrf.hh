#ifndef __LINALG_LAPACK_GEQRF_HH__
#define __LINALG_LAPACK_GEQRF_HH__

#include "matrix.hh"
#include "vector.hh"
#include "blas/dot.hh"
#include "blas/axpy.hh"
#include "operators.hh"

#include "openmp.hh"


namespace Linalg {
namespace Lapack {


/**
 * Internal function to calculate the matrix-vector product \f$a' = (I-2vv^T)a\f$.
 *
 * @ingroup lapack_internal
 */
inline void
__prod_householder(const Vector<double> &v, Vector<double> &a)
throw (ShapeError)
{
  double scale = -2*Blas::dot(v,a);
  Blas::__axpy(scale, v, a);
}


/**
 * Calculates the matrix-matrix product \f$A' = (I-2vv^T)A\f$.
 *
 * @ingroup lapack_internal
 */
inline void
__prod_householder(const Vector<double> &v, Matrix<double> &A)
throw (ShapeError)
{
  Vector<double> a;
  for(size_t j=0; j<A.cols(); j++) {
    a = A.col(j);
    __prod_householder(v, a);
  }
}



/**
 * Performs a QR decomposition with the same interface as LAPACKs geqrf function. In contrast to
 * geqrf, this variant uses Householder projectors and performs the decomposition on matrices in
 * any order (C or Fortran) while geqrf requires Fortran order.
 *
 * @param A Specifies a n-by-m matrix to be decomposed, with m <= n.
 * @param tau On exit, the upper triangular part of A will hold R and Q is encodes as elementary
 *        reflectors in the stricly lower-triangular part of A and in tau. dim(tau) >= m.
 * @throws ShapeError If m > n, n==0 or dim(tau) < m;
 *
 * @ingroup lapack
 */
void geqrf(Matrix<double> &A, Vector<double> &tau) throw (ShapeError)
{
  LINALG_SHAPE_ASSERT(A.rows()  >  0);
  LINALG_SHAPE_ASSERT(A.cols()  <= A.rows());
  LINALG_SHAPE_ASSERT(tau.dim() >= A.cols());

  size_t M = A.rows();
  size_t N = A.cols();

  // Iterate over all columns of A
  for (size_t i=0; i<N; i++) {
    // Create view on sub-matrix of A as Asub = A[i:,i:]
    // and sub-vector vsub = A[i:,i] = Asub[:,0];
    Matrix<double> Asub = A.sub(i,i, M-i,N-i);
    Vector<double> vsub = Asub.col(0);

    // Now, calculate Householder projector H = 1-2*vsub*vsub^T
    double alpha = vsub(0) < 0 ? std::abs(vsub) : -std::abs(vsub);
    vsub(0) += alpha; vsub /= std::abs(vsub);

    // Now, project all remaining columns of A[i:,i+1:], if some columns left
    for (size_t j=1; j<Asub.cols(); j++) {
      Vector<double> tmp = Asub.col(j);
      __prod_householder(vsub, tmp);
    }
    // Store v(0) in tau(i), alpha in A(i,i) and v[1:] in A[i+1:,i]
    tau(i) = vsub(0); vsub(0) = -alpha;
  }
}


#ifdef LINALG_HAS_OPENMP
/**
 * This function is identical to @c Lapack::geqrf, in constrast to that function, this function
 * uses ad-hoc parallelism by OpenMP.
 *
 * @param A Specifies a n-by-m matrix to be decomposed, with m <= n.
 * @param tau On exit, the upper triangular part of A will hold R and Q is encodes as elementary
 *        reflectors in the stricly lower-triangular part of A and in tau. dim(tau) >= m.
 * @param num_threads Specifies the number of threads to use. By default @c OpenMP::getMaxThreads()
 *        is used.
 * @throws ShapeError If m > n, n==0 or dim(tau) < m;
 *
 * @ingroup lapack
 */
void p_geqrf(Matrix<double> &A, Vector<double> &tau, size_t num_theads=OpenMP::getMaxThreads())
throw (ShapeError)
{
  LINALG_SHAPE_ASSERT(A.rows()  >  0);
  LINALG_SHAPE_ASSERT(A.cols()  <= A.rows());
  LINALG_SHAPE_ASSERT(tau.dim() >= A.cols());

  size_t M = A.rows();
  size_t N = A.cols();

  // Iterate over all columns of A
  for (size_t i=0; i<N; i++) {
    // Create view on sub-matrix of A as Asub = A[i:,i:]
    // and sub-vector vsub = A[i:,i] = Asub[:,0];
    Matrix<double> Asub = A.sub(i,i, M-i,N-i);
    Vector<double> vsub = Asub.col(0);

    // Now, calculate Householder projector H = 1-2*vsub*vsub^T
    double alpha = vsub(0) < 0 ? std::abs(vsub) : -std::abs(vsub);
    vsub(0) += alpha; vsub /= std::abs(vsub);

    // Now, project all remaining columns of A[i:,i+1:], if some columns left
#pragma omp parallel for
    for (size_t j=1; j<Asub.cols(); j++) {
      Vector<double> tmp = Asub.col(j);
      __prod_householder(vsub, tmp);
    }
    // Store v(0) in tau(i), alpha in A(i,i) and v[1:] in A[i+1:,i]
    tau(i) = vsub(0); vsub(0) = -alpha;
  }
}
#endif

}
}
#endif // __LINALG_LAPACK_GEQRF_HH__
