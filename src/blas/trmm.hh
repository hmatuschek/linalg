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
  if (Bcol.transposed()) {
    Acol = Acol.t();
    Bcol = Bcol.t();
    left = !left;
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
  if ((Acol.isUpper() && Acol.transposed()) || (!Acol.isUpper() && ! Acol.transposed()))
    UPLO = 'L';
  if (Acol.transposed())
    TRANSA = 'T';
  if (Acol.isUnit())
    DIAG = 'U';

  // Call fortran function:
  dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, *Acol, &LDA, *Bcol, &LDB);
}

}
}

#endif // __LINALG_BLAS_TRMM_HH__
