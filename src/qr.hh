#ifndef __LINALG_QR_HH__
#define __LINALG_QR_HH__

#include <vector.hh>
#include <matrix.hh>


namespace Linalg {


/**
 * Performs a QR decomposition with the same interface as LAPACKs geqrf function. In contrast to
 * geqrf, this variant uses Householder projectors and performs the decomposition on matrices in
 * any order (C or Fortran).
 *
 * @todo This method is untested.
 */
inline void geqrf(Matrix<double> &A, Vector<double> &tau, Vector<double> &v)
{
  LINALG_SHAPE_ASSERT(A.rows() > 0);
  LINALG_SHAPE_ASSERT(A.cols() <= A.rows());
  LINALG_SHAPE_ASSERT(tau.dim() <= A.cols());

  size_t m = A.rows();
  size_t n = A.cols();

  for (size_t i=0; i<n-1; i++) {
    Matrix<double> Asub = A.sub(i,i, m-i,n-i);
    // Copy values if the i-th column as A[i,i:] into v
    v.sub(0, m-i).values() = Asub.col(i);
    // Now, calculate Householder projector as 1-tmp*tmp^T
    double alpha = v(0) < 0 ? -abs(v.sub(0, m-i)) : abs(v.sub(0, m-i));
    v(0) += alpha; v /= abs(v.sub(0, m-i));

    // Now, project all columns of A[:,i+1:]
    for (size_t j=0; j<n; j++) {

    }
  }
}

}

#endif // __LINALG_QR_HH__
