/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_LAPACK_TRTRI_HH__
#define __LINALG_LAPACK_TRTRI_HH__


/* Interface to Fortran function. */
extern "C" {
void dtrtri_(const char *UPLO, const char *DIAG, const int *N,
             double *A, const int *LDA,
             int *INFO);
}


#include "blas/utils.hh"
#include "trimatrix.hh"

namespace Linalg {
namespace Lapack {


/**
 * Interface to LAPACKs DTRTRI function.
 *
 * Computes in-place the inverse of the triangular matrix A.
 *
 * @throws Linalg::SingularMatrixError If one of the diagonal elements of A is zero.
 * @throws Linalg::LapackError This should not happen. Thrown if one of the arguments
 *         to LAPACKs DTRTRI is invalid.
 *
 * @ingroup lapack
 */
inline void trtri(TriMatrix<double> A) throw (SingularMatrixError, LapackError)
{
  // A must be quadratic:
  LINALG_SHAPE_ASSERT(A.rows() == A.cols());

  // Ensure, Acol is in column-major:
  char transa='N';
  TriMatrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, transa);

  // Get all the flags:
  char uplo = BLAS_UPLO_FLAG(Acol);
  char diag = BLAS_UNIT_DIAG_FLAG(Acol);
  int  N    = BLAS_NUM_ROWS(Acol, transa);
  int  lda  = BLAS_LEADING_DIMENSION(Acol);
  int  info = 0;

  // Call function
  dtrtri_(&uplo, &diag, &N, Acol.ptr(), &lda, &info);

  // Check of errors:
  if (0 == info)
    return;

  // Check for error:
  if (0 > info) {
    LapackError err;
    err << -info <<"-th argument to DTRTRI() has illegal value.";
    throw err;
  } else if (0 < info) {
    SingularMatrixError err;
    err << "Signular matrix: " << info << "-th diagonal element of A is 0.";
  }
}


}
}
#endif // TRTRI_HH
