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
#include "trimatrix.hh"
#include "symmatrix.hh"


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
inline void dpotrf(Matrix<double> &A, bool upper)
throw (ShapeError, IndefiniteMatrixError, LapackError)
{
  // Assert that A is quadratic:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());

  // Ensure column major:
  char uplo = upper ? 'U' : 'L';
  Matrix<double> Acol(A);
  if (1 == Acol.strides()[0]) {
    uplo = uplo == 'U' ? 'L' : 'U';
    Acol = Acol.t();
  }

  // Get flags:
  int  N    = Acol.rows();
  int  lda  = BLAS_LEADING_DIMENSION(Acol);
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


}
}
#endif // __LINALG_LAPACK_POTRF_HH__
