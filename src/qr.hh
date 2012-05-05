#ifndef __LINALG_QR_HH__
#define __LINALG_QR_HH__

#include <vector.hh>
#include <matrix.hh>


namespace Linalg {


/**
 * Calculates the matrix-vector product \f$a' = (I-2vv^T)a\f$.
 */
inline void __prod_householder(const Vector<double> &v, Vector<double> &a) throw (ShapeError)
{
  double scale = 2*dot(v,a);
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

    double scale = 2*dot(v,a);
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
  size_t T = std::min(M-1, N);

  for (size_t i=0; i<T; i++) {
    // Create view on sub-matrix of A as Asub = A[i:,i:]
    // and sub-vector of v as vsub = v[:-i]
    Matrix<double> Asub = A.sub(i,i, M-i,N-i);
    Vector<double> vsub = v.sub(0,M-i);

    // Copy values if the i-th column as A[i:,i] into v
    vsub.values() = Asub.col(i);
    // Now, calculate Householder projector as 1-tmp*tmp^T
    double alpha = v(0) < 0 ? -abs(vsub) : abs(vsub);
    vsub(0) += alpha; vsub /= abs(vsub);

    // Store v(0) in tau(i), alpha in A(i,i) and v[1:] in A[i+1:,i]
    tau(i) = v(0); Asub.col(0).values() = vsub; Asub(0,0) = alpha;
    // Now, project all remaining columns of A[i:,i+1:], if some columns left
    for (size_t j=i+1; j<N; j++) {
      __prod_householder(vsub, Asub.col(j));
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
void ormqr(const Matrix<double> &A, Vector<double> &tau, Matrix<double> &B,
           bool trans=false, bool left=true) throw (ShapeError)
{
  // For all columns in B:
  for (size_t i=0; i<B.cols(); i++)
    ormqr(A, tau, B.col(i), trans, left);
}


/**
 * \todo This function is untested.
 */
void ormqr(const Matrix<double> &A, Vector<double> &tau, Vector<double> &b,
           bool trans=false, bool left=true) throw (ShapeError)
{

}


}

#endif // __LINALG_QR_HH__
