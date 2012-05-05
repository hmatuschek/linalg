#ifndef __LINALG_QR_HH__
#define __LINALG_QR_HH__

#include "vector.hh"
#include "matrix.hh"
#include "blas/dot.hh"
#include "operators.hh"



namespace Linalg {


/**
 * Internal function to calculate the matrix-vector product \f$a' = (I-2vv^T)a\f$.
 */
inline void __prod_householder(const Vector<double> &v, Vector<double> &a) throw (ShapeError)
{
  double scale = 2*Blas::dot(v,a);
  for (size_t i=0; i<v.dim(); i++)
    a(i) -= scale*v(i);
}


/**
 * Calculates the matrix-matrix product \f$A' = (I-2vv^T)A\f$.
 */
inline void __prod_householder(const Vector<double> &v, Matrix<double> &A) throw (ShapeError)
{
  for(size_t j=0; j<A.cols(); j++) {
    Vector<double> a = A.col(j);

    double scale = 2*Blas::dot(v,a);
    for (size_t i=0; i<v.dim(); i++)
      a(i) -= scale*v(i);
  }
}



/**
 * Performs a QR decomposition with the same interface as LAPACKs geqrf function. In contrast to
 * geqrf, this variant uses Householder projectors and performs the decomposition on matrices in
 * any order (C or Fortran) while geqrf requires Fortran order.
 *
 * @todo This method is untested.
 */
void geqrf(Matrix<double> &A, Vector<double> &tau, Vector<double> &v) throw (ShapeError)
{
  LINALG_SHAPE_ASSERT(A.rows()  >  0);
  LINALG_SHAPE_ASSERT(A.cols()  <= A.rows());
  LINALG_SHAPE_ASSERT(tau.dim() >= A.cols());
  LINALG_SHAPE_ASSERT(v.dim()   >= A.rows());

  size_t M = A.rows();
  size_t N = A.cols();

  for (size_t i=0; i<N; i++) {
    // Create view on sub-matrix of A as Asub = A[i:,i:]
    // and sub-vector of v as vsub = v[:-i]
    Matrix<double> Asub = A.sub(i,i, M-i,N-i);
    Vector<double> vsub = v.sub(0,M-i);

    // Copy values if the i-th column as A[i:,i] into v
    vsub.values() = Asub.col(i);
    // Now, calculate Householder projector as 1-tmp*tmp^T
    double alpha = vsub(0) < 0 ? -std::abs(vsub) : std::abs(vsub);
    vsub(0) += alpha; vsub /= std::abs(vsub);

    // Store v(0) in tau(i), alpha in A(i,i) and v[1:] in A[i+1:,i]
    tau(i) = vsub(0); Asub.col(0).values() = vsub; Asub(0,0) = alpha;
    // Now, project all remaining columns of A[i:,i+1:], if some columns left
    for (size_t j=1; j<Asub.cols(); j++) {
      Vector<double> tmp = Asub.col(j);
      __prod_householder(vsub, tmp);
    }
  }
}


/**
 * Simpler interface to @c geqrf.
 */
void geqrf(Matrix<double> &A, Vector<double> &tau) throw (ShapeError)
{
  // allocate some working memory...
  Vector<double> v(A.rows());
  geqrf(A, tau, v);
}



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

#endif // __LINALG_QR_HH__
