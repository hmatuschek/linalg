#ifndef __LINALG_LAPACK_GETC2_HH__
#define __LINALG_LAPACK_GETC2_HH__

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
 * Nice interface to the LAPACK function "dgetc2". Computes the pivoting for a LU decomposition.
 *
 * @param a Specifies a N-by-N real matrix (in/out).
 * @param ipiv N-dim Vector<int> of indices of rows (out).
 * @param jpiv N-dim Vector<int> of indices of columns (out).
 */
void getc2(Matrix<double> &a, Vector<int> &ipiv, Vector<int> &jpiv, size_t &rank)
{
  // Assert A being a N-by-N matrix:
  if (a.rows() != a.cols())
  {
    Linalg::IndexError err;
    err << "a must be a quadratic matrix!";
    throw err;
  }

  int n   = a.rows();
  int lda = a.leading_dimension();
  int info = 0;

  // Calling fortran:
  dgetc2_(&n, *a, &lda, *ipiv, *jpiv, &info);

  // Check for errors:
  if (0 > info)
  {
    Linalg::LapackError err("Lapack: ");
    err << "LAPACK returned error: " << -info;
    throw err;
  }

  rank = n-info;
}


}
}

#endif
