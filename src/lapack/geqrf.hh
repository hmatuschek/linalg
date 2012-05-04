#ifndef __LINALG_LAPACK_GEQRF_HH__
#define __LINALG_LAPACK_GEQRF_HH__

#include "matrix.hh"
#include "workspace.hh"


namespace Linalg {
namespace Lapack {

/**
 * Interface to LAPACKs DGEQRF() function, computing the QR-decomposition of a general matrix A.
 *
 * @param[in,out] A Specifies the matrix to compute the QR decompotion.
 * @param[out] tau Will hold the factors of the elementary reflectors. Must be at least of length
 *                 @c min(A.rows(), A.cols()).
 * @param[in] work Specifies the work-space instance to use.
 *
 * @todo untested.
 *
 * @ingroup lapack
 */
void geqrf(Matrix<double> &A, Vector<double> &tau, Workspace &work) throw (ShapeError, LapackError)
{
  // Check if length of tau >= min(cols(A), rows(A))
  LINALG_SHAPE_ASSERT(tau.dim() >= std::min(A.rows(), A.cols()));

  // Ensure column-major form:
  char transa = 'N';
  Matrix<double> Acol = A; BLAS_ENSURE_COLUMN_MAJOR(Acol, transa);

  // Allocate some work-space
  double *work_ptr = work.ensure<double>(A.cols());
  int lwork = work.size<double>();

  int M   = BLAS_NUM_ROWS(Acol);
  int N   = BLAS_NUM_COLS(Acol);
  int lda = BLAS_LEADING_DIMENSION(Acol);
  int info = 0;

  dgeqrf_(&M, &N, Acol.ptr(), &lda, tau.ptr(), work_ptr, &lwork, &info);

  if (0 > info) {
    LapackError err;
    err << "In geqrf(): " << -i << "-th argument to DGEQRF() had illegal value.";
    throw err;
  }
}


}
}
#endif // __LINALG_LAPACK_GEQRF_HH__
