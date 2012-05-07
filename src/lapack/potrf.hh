/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_LAPACK_POTRF_HH__
#define __LINALG_LAPACK_POTRF_HH__


/* Interface to Fortran function. */
extern "C" {
void dpotrf_(const char *UPLO, const int *N, double *A, const int *LDA, int *INFO);
}


#include "blas/utils.hh"
#include "utils.hh"


namespace Linalg {
namespace Lapack {


/**
 * Interface to LAPACKs DPOTRF function, computing the Cholesky decomposition of a
 * real symmetric matrix.
 *
 * @param A Specifies the symetric matrix, only the upper or lower part is used.
 * @param upper If true, the upper-triangular part of the symmetric matrix is stored in A.
 *
 * @ingroup lapack
 */
inline void
dpotrf(Matrix<double> &A, bool upper)
throw (ShapeError, IndefiniteMatrixError, LapackError)
{
  // Assert that A is quadratic:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());

  // Ensure column major:
  char uplo = upper ? 'U' : 'L';
  Matrix<double> Acol(A);
  if (1 == Acol.strides()[1]) {
    uplo = (uplo == 'U') ? 'L' : 'U';
    Acol = Acol.t();
  }

  // Get flags:
  int  N    = Acol.rows();
  int  lda  = Acol.strides()[1];
  int  info = 0;

  // Call function
  dpotrf_(&uplo, &N, Acol.ptr(), &lda, &info);

  // Check for errors:
  if (0 != info) {
    if (0 > info) {
      LapackError err;
      err << "Argument error: " << -info << "-th argument to DPOTRF() has illegal value.";
      throw err;
    } else {
      IndefiniteMatrixError err;
      err << "The leading minor of order " << info
          << " is not positive definite, can not compute Cholesky decomposition.";
      throw err;
    }
  }
}



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
 * This algorithm is identical to the one implemented in @c __potrf_banachiewicz but in contrast
 * to that algorithm, this one operates on the lower-triangualar.
 *
 * @ingroup lapack_internal
 */
template <class Scalar>
void
__potrf_crout(Matrix<Scalar> &A)
throw (IndefiniteMatrixError)
{
  size_t i,j;
  for (j=0; j<A.cols(); j++)
  {
    for (size_t k=0; k<j; k++)
      A(j,j) -= __prod_aacc(A(j,k));

    // Check if the j-th diagonal is positve and real:
    if (! __is_real_pos(A(j,j)) )
    {
      IndefiniteMatrixError err; err << "Can not compute Cholesky dec. of A, "
                                     << j << "-th diagonal element is <= 0!";
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
 * This algorithm is identical to the one implemented in @c __potrf_crout but in contrast
 * to that algorithm, this one operates on the upper-triangualar.
 *
 * @ingroup lapack_internal
 */
template <class Scalar>
void
__potrf_banachiewicz(Matrix<Scalar> &A)
throw (IndefiniteMatrixError)
{
  size_t i,j, N=A.cols();
  for (i=0; i<N; i++)
  {
    for (size_t k=0; k<i; k++)
      A(i,i) -= __prod_aacc(A(k,i));

    // Check if the j-th diagonal is positve and real:
    if (! __is_real_pos(A(i,i)) )
    {
      IndefiniteMatrixError err; err << "Can not compute Cholesky dec. of A, "
                                     << i << "-th diagonal element is <= 0!";
      throw err;
    }
    A(i,i) = std::sqrt(__get_real(A(i,i)));

    // calculate off-diagonal elements
    for (j=i+1; j<N; j++) {
      for (size_t k=0; k<i; k++)
        A(i,j) -= __prod_abcc(A(k,j), A(k,i));

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
 *
 * @throws IndefiniteMatrixError If one of the diagonal elements are <= 0.
 *
 * @ingroup lapack
 */
template <class Scalar>
inline void
potrf(Matrix<Scalar> &A, bool upper)
throw (IndefiniteMatrixError)
{
  // Check if A is square:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());

  // Dispatch
  if (upper) {
    __potrf_banachiewicz(A);
  } else {
    __potrf_crout(A);
  }
}


}
}
#endif // __LINALG_LAPACK_POTRF_HH__
