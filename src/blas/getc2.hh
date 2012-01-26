/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_GETC2_HH__
#define __LINALG_BLAS_GETC2_HH__

#include "blas/utils.hh"
#include "matrix.hh"
#include "vector.hh"

/*
 * Interface to Fortran function.
 */
extern "C" {
extern void dgetc2_(const int *n, double *a, const int *lda,
                    int *ipiv, int *jpiv, int *info);
}



namespace Linalg {
namespace Lapack {


/**
 * Nice interface to the BLAS function "dgetc2". Computes the full pivoting LU-decomposition:
 *
 * \f[ A = P \cdot L \cdot U \cdot Q \f]
 *
 * @todo This function is untested!
 *
 * @param a Specifies a N-by-N real matrix (in/out).
 * @param ipiv N-dim Vector<int> of indices of rows (out).
 * @param jpiv N-dim Vector<int> of indices of columns (out).
 *
 * @ingroup blas2
 */
void getc2(Matrix<double> &A, Vector<int> &ipiv, Vector<int> &jpiv, size_t &rank)
{
  // Get A in column-major form:
  Matrix<double> Acol = BLAS_TO_COLUMN_MAJOR(A);

  // Assert A being a N-by-N matrix:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());
  LINALG_SHAPE_ASSERT(Acol.rows() == ipiv.dim());
  LINALG_SHAPE_ASSERT(Acol.rows() == jpiv.dim());

  int N   = BLAS_NUM_ROWS(Acol);
  int lda = BLAS_LEADING_DIMENSION(Acol);
  int info = 0;

  int *ipiv_ptr = *ipiv;
  int *jpiv_ptr = *jpiv;

  if (BLAS_IS_TRANSPOSED(Acol)) {
    Acol = Acol.t();
    std::swap(ipiv_ptr, jpiv_ptr);
  }

  // Calling fortran:
  dgetc2_(&N, *A, &lda, ipiv_ptr, jpiv_ptr, &info);

  // Check for errors:
  if (0 > info)
  {
    Linalg::LapackError err("Lapack: ");
    err << "DGETC2 returned error: " << -info;
    throw err;
  }

  rank = N-info;
}


}
}

#endif
