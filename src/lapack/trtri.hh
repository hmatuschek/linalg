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
 * @todo This function is untested yet.
 *
 * @ingroup lapack
 */
inline void trtri(TriMatrix<double> &A)
{
  TriMatrix<double> Acol = BLAS_TO_COLUMN_MAJOR(A);

  // A must be quadratic:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());

  // Get all the flags:
  char UPLO = BLAS_UPLO_FLAG(Acol);
  char DIAG = BLAS_UNIT_DIAG_FLAG(Acol);
  int  N    = BLAS_NUM_ROWS(Acol);
  int  LDA  = BLAS_LEADING_DIMENSION(Acol);
  int  INFO = 0;

  // Call function
  dtrtri_(&UPLO, &DIAG, &N, *Acol, &LDA, &INFO);

  // Check of errors:
  if (0 == INFO)
    return;

  if (0 > INFO) {
    RuntimeError err;
    err << -INFO <<"-th argument to DTRTRI() has illegal value.";
    throw err;
  } else if (0 < INFO) {
    SingularMatrixError err;
    err << "Signular matrix: " << INFO << "-th diagonal element of A is 0.";
  }
}


}
}
#endif // TRTRI_HH
