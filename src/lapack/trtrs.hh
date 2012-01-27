#ifndef __LINALG_LAPACK_TRTRS_HH__
#define __LINALG_LAPACK_TRTRS_HH__

#include "blas/utils.hh"
#include "blas/trsm.hh"
#include "matrix.hh"
#include "trimatrix.hh"


namespace Linalg {
namespace Lapack {

/**
 * Reimplementation of LAPACKs ?TRTRS.
 *
 * Solves the triangular system:
 * \f[A\cdot X = B\f]
 *
 * The LAPACK ?TRTRS implementation was not flexible enougth, so I decided to reimplement
 * this function as a template. The LAPACK function simply checks if there is a 0 on the
 * diagonal of @c A. If not it calls BLASs ?TRSM, so do I.
 *
 * See also: @c Blas::trsm
 *
 * @ingroup lapack
 */
template <class Scalar>
void trtrs(const TriMatrix<Scalar> &A, Matrix<Scalar> &B)
{
  // Ensure column-major representation of matrices:
  TriMatrix<Scalar> Acol = BLAS_TO_COLUMN_MAJOR(A);
  Matrix<Scalar> Bcol = BLAS_TO_COLUMN_MAJOR(B);

  // A must be quadratic:
  LINALG_SHAPE_ASSERT(Acol.rows() == Acol.cols());
  LINALG_SHAPE_ASSERT(Acol.cols() == Bcol.rows());

  // If B is transposed, transpose A & B and swap sides:
  if (Bcol.isTransposed()) {
    Bcol = Bcol.t();
    Acol = Acol.t();
    left = !left;
  }

  // Check if A has zero on diagonal:
  if (! A.hasUnitDiag()) {
    for (size_t i=0; i<Acol.rows(); i++) {
      if (0 == Acol(i,i)) {
        SingulrMatrixError err;
        err << "Signular matrix: " << i << "-th diagonal element of A is 0!";
        throw err;
      }
    }
  }

  // Just call trsm() ...
  Blas::trsm(Acol, 1.0, Bcol, left);
}


}
}
#endif // TRTRS_HH
