/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_TRMM_HH__
#define __LINALG_BLAS_TRMM_HH__

// Iternface to Fortran function:
extern "C" {
void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double *a, int *lda, double *b, int *ldb);
}



#include "trimatrix.hh"


namespace Linalg {
namespace Blas {

/**
 * Implements an interface to the DTRMM function of BLAS.
 *
 * Calculates
 * \f[B = A\cdot B\f]
 * if @c left is true and
 * \f[B = B\cdot A\f]
 * if @c left is false.
 *
 * @todo This function is not tested yet.
 *
 * @param left Specifies if @c A is multiplied from left or right to @c B.
 * @param A Specifies a unit or non-unit, upper- or lower-triangular matrix.
 * @param B Specifies some general matrix.
 *
 * @ingroup blas3
 */
void trmm(bool left, double alpha, const TriMatrix<double> &A, Matrix<double> &B)
{
  // Make sure, A & B are in column-major form:
  TriMatrix<double> Acol = A.colMajor();
  Matrix<double> Bcol = B.colMajor();

  // If Bcol is transposed: swap side and transpose Acol:
  if (Bcol.isTransposed()) {
    Acol = Acol.t();
    Bcol = Bcol.t();
    left = !left;
  }

  // Check dimensions:
  if (left) {
    LINALG_SHAPE_ASSERT(Acol.rows() == Bcol.rows());
    LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.rows());
  } else {
    LINALG_SHAPE_ASSERT(Acol.rows() == Bcol.cols());
    LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.cols());
  }

  char SIDE = 'L';
  char UPLO = 'U';
  char TRANSA = 'N';
  char DIAG = 'N';
  int  M = Bcol.rows();
  int  N = Bcol.cols();
  double ALPHA = alpha;
  int LDA = A.stride();
  int LDB = B.stride();

  if (!left)
    SIDE = 'R';
  if ((!Acol.isUpper() && !Acol.isTransposed()) || (Acol.isUpper() && Acol.isTransposed()))
    UPLO = 'L';
  if (Acol.isTransposed())
    TRANSA = 'T';
  if (Acol.isUnit())
    DIAG = 'U';

  // Call fortran function:
  dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, *Acol, &LDA, *Bcol, &LDB);
}

}
}

#endif // __LINALG_BLAS_TRMM_HH__
