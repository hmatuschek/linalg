#ifndef __LINALG_BLAS_TRSM_HH__
#define __LINALG_BLAS_TRSM_HH__


extern "C" {
void trsm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG,
           const int *M, const int *N,
           const double *alpha, const double *A, const int *LDA,
           double *B, const int *LDB);
}


#include "blas/utils.hh"
#include "matrix.hh"
#include "trimatrix.hh"


namespace Linalg {
namespace Blas{

/**
 * Interface to BLAS level 3 DTRSM function.
 *
 * Solves the triangular system \f$ A\cdot X = \alpha B\f$ if @c left=true
 * or \f$X \cdot A = \alpha B\f$ if @c left=false.
 *
 * @todo This function is untested yet.
 *
 * @ingroup blas3
 */
void trsm(const TriMatrix<double> &A, double alpha, Matrix<double> &B, bool left=true)
{
  TriMatrix<double> Acol = BLAS_TO_COLUMN_MAJOR(A);
  Matrix<double>    Bcol = BLAS_TO_COLUMN_MAJOR(B);

  if (Bcol.isTransposed()) {
    Acol = Acol.t();
    Bcol = Bcol.t();
    left = !left;
  }

  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());
  if (left) {
    LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.rows())
  } else {
    LINALG_SHAPE_ASSERT(Bcol.cols() == Acol.rows());
  }

  char SIDE   = left ? 'L' : 'R';
  char UPLO   = BLAS_UPLO_FLAG(Acol);
  char TRANSA = BLAS_TRANSPOSED_FLAG(Acol);
  char DIAG   = BLAS_UNIT_DIAG_FLAG(Acol);
  int  M      = BLAS_NUM_ROWS(Bcol);
  int  N      = BLAS_NUM_COLS(Bcol);
  int  LDA    = BLAS_LEADING_DIMENSION(Acol);
  int  LDB    = BLAS_LEADING_DIMENSION(Bcol);

  dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &alpha, *Acol, &LDA, *Bcol, &LDB);
}

}
}
#endif // TRSM_HH
