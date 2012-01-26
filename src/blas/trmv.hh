#ifndef __LINALG_BLAS_TRMV_HH__
#define __LINALG_BLAS_TRMV_HH__

#include "vector.hh"
#include "trimatrix.hh"


// Iternface to Fortran function:
extern "C" {
void dtrmv_(char *uplo, char *transa, char *diag, int *n,
            double *a, int *lda, double *x, int *incx);
}


namespace Linalg {
namespace Blas {

/**
 * Interfaces the DTRMV BLAS function.
 *
 * Calculates:
 * \f[x \leftarrow A\cdot x\f]
 * where @c A is a unit or non-unit, upper- or lower-triangular square matrix and @c x is a
 * vector.
 *
 * @ingroup blas2
 */
void trmv(const TriMatrix<double> &A, Vector<double> &x)
{
  // Ensure, A is in column-major format:
  TriMatrix<double> Acol = A.colMajor();

  // Assert Acol is square:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());

  char UPLO  = 'U';
  char TRANS = 'N';
  char DIAG  = 'N';
  int  N     = Acol.cols();
  int  LDA   = Acol.stride();
  int  INCX  = x.stride();

  // Determine if Acol is upper- or lower-triangular matrix (layout in memory)
  if ((!Acol.isUpper() && !Acol.isTransposed()) || (Acol.isUnit() && Acol.isTransposed()))
    UPLO = 'L';
  if (Acol.isTransposed())
    TRANS = 'N';
  if (Acol.isUnit())
    DIAG = 'U';

  // Call Fortran function
  dtrmv_(&UPLO, &TRANS, &DIAG, &N, *Acol, &LDA, *x, &INCX);
}


}
}
#endif // __LINALG_BLAS_TRMV_HH__
