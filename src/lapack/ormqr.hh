#ifndef __LINALG_LAPACK_ORMQR_HH__
#define __LINALG_LAPACK_ORMQR_HH__

#include "geqrf.hh"


namespace Linalg {
namespace Lapack {


/**
 * \todo This function is untested.
 */
void ormqr(const Matrix<double> &A, const Vector<double> &tau, Vector<double> &b,
           Vector<double> &v, bool trans=false, bool left=true) throw (ShapeError)
{
  LINALG_SHAPE_ASSERT(A.rows() == b.dim());
  LINALG_SHAPE_ASSERT(v.dim() >= A.rows());

  size_t M = A.rows();
  size_t N = A.cols();
  size_t K = std::min(M-1, N);

  // Apply Householder projectors in reverse order
  if ( (left && !trans) || (!left && trans)) {
    for (int i=K; i>=0; i--) {
      // Assemble Householder projector in v:
      v.sub(1,M-i-1).values() = A.col(i).sub(i+1, M-i-1); v(0) = tau(i);
      // Apply Householder
      Vector<double> b_sub = b.sub(i, M-i);
      __prod_householder(v.sub(0, M-i), b_sub);
    }
  }
  // Apply Householder projectors in normal order
  else {
    for (size_t i=0; i<K; i++) {
      // Assemble Householder projector in v:
      v.sub(1,M-i-1).values() = A.col(i).sub(i+1, M-i-1); v(0) = tau(i);
      // Apply Householder
      Vector<double> b_sub(b.sub(i, M-i));
      __prod_householder(v.sub(0, M-i), b_sub);
    }
  }
}


/**
 * \todo This function is untested.
 */
void ormqr(const Matrix<double> &A, const Vector<double> &tau, Matrix<double> &B,
           Vector<double> &v, bool trans=false, bool left=true) throw (ShapeError)
{
  LINALG_SHAPE_ASSERT(A.rows() > 0);
  LINALG_SHAPE_ASSERT(A.rows() >= A.cols());
  LINALG_SHAPE_ASSERT(A.cols() <= tau.dim());

  // For all columns in B:
  for (size_t i=0; i<B.cols(); i++) {
    Vector<double> Bcol;
    if ((left && !trans) || (!left && trans))
      Bcol = B.col(i);
    else
      Bcol = B.row(i);
    ormqr(A, tau, Bcol, v, trans, left);
  }
}


}
}
#endif // __LINALG_LAPACK_ORMQR_HH__
