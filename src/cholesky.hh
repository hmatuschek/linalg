#ifndef __LINALG_CHOLESKY_HH__
#define __LINALG_CHOLESKY_HH__

#include <matrix.hh>
#include <utils.hh>


namespace Linalg {

/**
 * This function implements the Cholesky-Crout algorithm, calculating in-place the Cholesky
 * decomposition of a real-symmetric or complex-hermitian matrix stored in the lower-triangular
 * part of A.
 *
 * \f[
 *  L_{jj} &=& \sqrt{A_{j,j}-\sum_{k=0}^{j-1}{L_{j,k}L_{j,k}^*}} \\
 *  L_{ij} &=& \frac{1}{L_{jj}}\left(A_{ij}-\sum_{k=1}^{j-1}L_{ik}L_{jk}^*\right)\quad, i>j
 * \f]
 *
 * This algorithm is identical to the one implemented in @c __cholesky_banachiewicz but in contrast
 * to that algorithm, this one operates on the lower-triangualar.
 */
template <class Scalar>
void __cholesky_crout(Matrix<Scalar> &A) throw (SingularMatrixError)
{
  size_t i,j;
  for (j=0; j<A.cols(); j++)
  {
    for (size_t k=0; k<j; k++)
      A(j,j) -= __prod_aacc(A(j,k));

    // Check if the j-th diagonal is positve and real:
    if (! __is_real_pos(A(j,j)) )
    {
      SingularMatrixError err; err << "Can not compute Cholesky dec. of A, "
                                   << j << "-th diagonal element is < 0!";
      throw err;
    }
    A(j,j) = std::sqrt(__get_real(A(j,j)));

    i = j+1;
    for (; i<A.rows(); i++) {
      for (size_t k=0; k<j; k++)
        A(i,j) -= __prod_abcc(A(i,k), A(j,k));
      A(i,j) /= A(j,j);
    }
  }
}



/**
 * This function implements the Cholesky-Banachiewicz algorithm, calculating in-place the Cholesky
 * decomposition of a real-symmetric or complex-hermitian matrix stored in the upper-triangular
 * part of A.
 *
 * This algorithm is identical to the one implemented in @c __cholesky_crout but in contrast
 * to that algorithm, this one operates on the upper-triangualar.
 */
template <class Scalar>
void __cholesky_banachiewicz(Matrix<Scalar> &A) throw (SingularMatrixError)
{
  size_t i,j, N=A.cols();
  for (i=0; i<N; i++)
  {
    for (size_t k=0; k<i; k++)
      A(i,i) -= __prod_aacc(A(k,i));

    // Check if the j-th diagonal is positve and real:
    if (! __is_real_pos(A(i,i)) )
    {
      SingularMatrixError err; err << "Can not compute Cholesky dec. of A, "
                                   << i << "-th diagonal element is < 0!";
      throw err;
    }
    A(i,i) = std::sqrt(__get_real(A(i,i)));

    // calculate off-diagonal elements
    for (j=i+1; j<N; j++) {
      for (size_t k=0; k<i; k++)
        A(i,j) -= __prod_abcc(A(k,i), A(k,j));

      A(i,j) /= A(i,i);
    }
  }
}



/**
 * This function calculates the Cholesky decomposition of a real-symmetric or complex-hermitan
 * matrix stored in the upper or lower triangular part of A. The result is stored into the
 * upper or lower triangular part of the matrix.
 *
 * @param A Holds the upper or lower part of the real-symmetric or complex-hermitic matrix.
 * @param upper If true, the upper triangular part of A is given.
 * @throws SingularMatrixError If one of the diagonal elements are 0.
 */
template <class Scalar>
inline void cholesky(Matrix<Scalar> &A, bool upper) throw (SingularMatrixError)
{
  // Check if A is square:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());

  // Dispatch
  if (upper) {
    __cholesky_banachiewicz(A);
  } else {
    __cholesky_crout(A);
  }
}

}

#endif // __LINALG_CHOLESKY_HH__
