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
 * real- symmetric matrix.
 *
 * @note On extit, the Cholesky decomposition is stored into A and a @c TriMatrix view is returned,
 *       referring to the upper- or lower-triangular part of A holding the decomposition.
 *
 * @ingroup lapack
 */
inline TriMatrix<double> potrf(SymMatrix<double> A)
throw (ShapeError, IndefiniteMatrixError, LapackError)
{
  // Assert that A is quadratic:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());

  // Ensure column major:
  char transa='N';
  SymMatrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, transa);

  // Get flags:
  char uplo = BLAS_UPLO_FLAG(Acol);
  int  N    = BLAS_NUM_ROWS(Acol, transa);
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

  // Assemble triangular matrix view
  return TriMatrix<double>(Acol, Acol.isUpper(), false);
}


}
}
#endif // __LINALG_LAPACK_POTRF_HH__
