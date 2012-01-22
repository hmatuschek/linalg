#ifndef __FLUC_LINALG_BLAS_HH__
#define __FLUC_LINALG_BLAS_HH__

#include "workspace.hh"
#include "arrayview.hh"


namespace Fluc {
namespace Linalg {


/**
 * Direct access to the dgemv LAPACK function.
 *
 * @ingroup linalg
 */
void gemv(const double alpha, const Matrix<double> &a, const Vector<double> &x,
          const double beta, Vector<double> &y);

/**
 * Direct access to the dgemm LAPACK function.
 *
 * @ingroup linalg
 */
void gemm(const double alpha, const Matrix<double> &a, const Matrix<double> &b,
          const double beta, Matrix<double> &c);


/**
 * Implements the matrix-matrix product.
 *
 * @ingroup linalg
 */
Matrix<double> dot(const Matrix<double> &a, const Matrix<double> &b);

/**
 * Implements the matrix-vector product.
 *
 * @ingroup linalg
 */
Vector<double> dot(const Matrix<double> &a, const Vector<double> &b);


}
}
#endif // __FLUC_LINALG_BLAS_HH__
