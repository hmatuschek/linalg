#ifndef __LINALG_BLAS_GETC2_HH__
#define __LINALG_BLAS_GETC2_HH__

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
 * @param a Specifies a N-by-N real matrix (in/out).
 * @param ipiv N-dim Vector<int> of indices of rows (out).
 * @param jpiv N-dim Vector<int> of indices of columns (out).
 *
 * @ingroup blas2
 */
void getc2(Matrix<double> &A, Vector<int> &ipiv, Vector<int> &jpiv, size_t &rank)
{
  // Get A in column-major form:
  Matrix<double> Acol = A.colMajor();

  // Assert A being a N-by-N matrix:
  if (Acol.rows() != Acol.cols())
  {
    Linalg::IndexError err;
    err << "A must be a quadratic matrix!";
    throw err;
  }

  if (Acol.rows() != ipiv.dim() || Acol.rows() != jpiv.dim())
  {
    Linalg::IndexError err;
    err << "rows(A) == cols(A) == dim(ipiv) == dim(jpiv) not satisfied!";
    throw err;
  }

  int N   = Acol.rows();
  int lda = Acol.stride();
  int info = 0;

  double *ipiv_ptr = *ipiv;
  double *jpiv_ptr = *jpiv;

  if (Acol.transposed())
    std::swap(ipiv_ptr, jpiv_ptr);

  // Calling fortran:
  dgetc2_(&N, *A, &lda, ipiv_ptr, jpiv_ptr, &info);

  // Check for errors:
  if (0 > info)
  {
    Linalg::LapackError err("Lapack: ");
    err << "DGETC2 returned error: " << -info;
    throw err;
  }

  rank = n-info;
}


}
}

#endif
