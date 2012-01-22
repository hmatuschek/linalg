#ifndef __FLUC_LINALG_LAPACK_HH__
#define __FLUC_LINALG_LAPACK_HH__

#include "arrayview.hh"


namespace Fluc {
namespace Linalg {


/**
 * Nice interface to the LAPACK function "dgetc2". Computes the pivoting for a LU decomposition.
 *
 * @param a Specifies a N-by-N real matrix (in/out).
 * @param ipiv N-dim Vector<int> of indices of rows (out).
 * @param jpiv N-dim Vector<int> of indices of columns (out).
 */
void getc2(Matrix<double> &a, Vector<int> &ipiv, Vector<int> &jpiv, size_t &rank);


/**
 * Computes the complete pivoted LU decomposition of the N-by-M matrix
 * \f$A = P\cdot L\cdot U\cdot Q\f$ where P,Q are row/column permutation matrices, L being a lower
 * unit-triangular matrix and U being a upper-triangular matrix.
 *
 * @ingroup linalg
 */
class LUDec
{
protected:
  size_t rank;

  /**
   * Holds the L and U matrices, the unit diagonal of L is not stored.
   */
  Matrix<double> lu;

  /**
   * Index mapping of rows i -> ipiv(i).
   */
  Vector<int> ipiv;

  /**
   * Index mapping of columns j -> jpiv(j).
   */
  Vector<int> jpiv;


public:
  /**
   * Performs the LU decomposition of a.
   */
  LUDec(const Matrix<double> &a);

  /**
   * Returns the permutation matrix P.
   */
  Matrix<double> getP();

  /**
   * Returns the lower unit-triangular matrix L.
   */
  Matrix<double> getL();

  /**
   * Returns the upper-triangular matrix U.
   */
  Matrix<double> getU();

  /**
   * Returns the permutation matrix Q.
   */
  Matrix<double> getQ();

  Vector<int> getRowPerm();
  Vector<int> getColPerm();
};


}
}

#endif // LAPACK_HH
