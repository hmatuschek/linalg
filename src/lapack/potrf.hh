#ifndef __LINALG_LAPACK_POTRF_HH__
#define __LINALG_LAPACK_POTRF_HH__


extern "C" {
void dportf_(const char *UPLO, const int *N, double *A, const int *LDA, int *INFO);
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
 *       referring to the upper- or lower-triangular decomposition.
 * @ingroup lapack
 */
TriMatrix<double> potrf(SymMatrix<double> &A)
{
  SymMatrix<double> Acol = BLAS_TO_COLUMN_MAJOR(A);

  // Assert that A is quadratic:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());

  char UPLO = BLAS_UPLO_FLAG(Acol);
  int  N    = BLAS_NUM_ROWS(Acol);
  int  LDA  = BLAS_LEADING_DIMENSION(Acol);
  int  INFO = 0;

  dpotrf_(&UPLO, &N, *Acol, &LDA, &INFO);

  if (0 == INFO)
    return;

  if (0 > INFO) {
    RuntimeError err;
    err << "Argument error: " << -INFO << "-th argument to DPOTRF() has illegal value.";
    throw err;
  } else {
    IndefiniteMatrixError err;
    err << "The leading minor of order " << INFO
        << " is not positive definite, can not compute Cholesky decomposition.";
    throw err;
  }

  return TriMatrix<double>(Acol, Acol.isUpper(), false);
}


}
}
#endif // __LINALG_LAPACK_POTRF_HH__
