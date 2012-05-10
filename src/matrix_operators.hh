#ifndef __LINALG_MATRIX_OPERATORS_HH__
#define __LINALG_MATRIX_OPERATORS_HH__

#include "matrix.hh"
#include "blas/gemm.hh"
#include "blas/gemv.hh"


namespace Linalg {

/**
 * Implements the common matrix-matrix product.
 *
 * @note This operator is quiet inefficient, as a matrix is allocated to hold the itermediate result
 *       of the operation. Whenever possible, use the in-place operator *=.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Matrix<Scalar>
operator* (const Matrix<Scalar> &lhs, const Matrix<Scalar> &rhs) throw (ShapeError)
{
  // Allocate matrix for result:
  Matrix<Scalar> result(lhs.rows(), rhs.cols());
  // Use Blas::gemm() to compute product
  Blas::gemm(Scalar(1), lhs, rhs, Scalar(0), result);
  // Pass ownership of result-matrix to caller...
  return result;
}



/**
 * Implements the common matrix-vector product.
 *
 * @note This operator is quiet inefficient, as a vector is allocated to hold the itermediate result
 *       of the operation. Whenever possible, use the in-place operator *=.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Vector<Scalar>
operator* (const Matrix<Scalar> &lhs, const Vector<Scalar> &rhs) throw (ShapeError)
{
  // Allocate vector for result:
  Vector<Scalar> result(rhs.dim());
  // Use Blas::gemv() to compute product
  Blas::gemv(Scalar(1), lhs, rhs, Scalar(0), result);
  // Pass ownership of result vector to caller...
  return result;
}


/**
 * Implements the common matrix sum.
 *
 * @note This operator is quiet inefficient, as a matrix is allocated to hold the itermediate result
 *       of the operation. Whenever possible, use the in-place operator +=.
 *
 * @ingroup operators
 */
inline Matrix<double>
operator+ (const Matrix<double> &lhs, const Matrix<double> &rhs)
{
  // Check shape:
  LINALG_SHAPE_ASSERT(lhs.rows() == rhs.rows());
  LINALG_SHAPE_ASSERT(rhs.cols() == rhs.cols());

  // Allocate result matrix.
  Matrix<double> result(lhs.rows(), lhs.cols());

  // Perform sum:
  for (size_t i=0; i<lhs.rows(); i++) {
    for (size_t j=0; j<lhs.cols(); j++) {
      result(i,j) = lhs(i,j) + rhs(i,j);
    }
  }

  // Pass ownership of result to caller:
  return result;
}


/**
 * Implements the common matrix difference.
 *
 * @note This operator is quiet inefficient, as a matrix is allocated to hold the itermediate result
 *       of the operation. Whenever possible, use the in-place operator -=.
 *
 * @ingroup operators
 */
inline Matrix<double>
operator- (const Matrix<double> &lhs, const Matrix<double> &rhs)
{
  // Check shape:
  LINALG_SHAPE_ASSERT(lhs.rows() == rhs.rows());
  LINALG_SHAPE_ASSERT(rhs.cols() == rhs.cols());

  // Allocate result matrix.
  Matrix<double> result(lhs.rows(), lhs.cols());

  // Perform sum:
  for (size_t i=0; i<lhs.rows(); i++) {
    for (size_t j=0; j<lhs.cols(); j++) {
      result(i,j) = lhs(i,j) - rhs(i,j);
    }
  }

  // Pass ownership of result to caller:
  return result;
}


/**
 * Implements the in-place matrix-scaleing.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Matrix<Scalar> & operator *= (Matrix<Scalar> &A, const Scalar &alpha) {
  // If column-major, work column by column:
  if (1 == A.strides()[0]) {
    for (size_t i=0; i<A.cols(); i++) {
      Vector<Scalar> col = A.col(i);
      Blas::scal(alpha, col);
    }
  } else {
    // Otherwise, row by row
    for (size_t i=0; i<A.rows(); i++) {
      Vector<Scalar> row = A.row(i);
      Blas::scal(alpha, row);
    }
  }

  return A;
}


/**
 * Implements the in-place matrix-scaleing.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Matrix<Scalar> & operator /= (Matrix<Scalar> &A, const Scalar &alpha) {
  // If column-major, work column by column:
  if (1 == A.strides()[0]) {
    for (size_t i=0; i<A.cols(); i++) {
      Vector<Scalar> col = A.col(i);
      Blas::scal(Scalar(1)/alpha, col);
    }
  } else {
    // Otherwise, row by row
    for (size_t i=0; i<A.rows(); i++) {
      Vector<Scalar> row = A.row(i);
      Blas::scal(Scalar(1)/alpha, row);
    }
  }

  return A;
}


/**
 * Output operator.
 *
 * @ingroup operators
 */
template <class Scalar>
std::ostream &operator<< (std::ostream &stream, Matrix<Scalar> &matrix)
{
  stream << "[";
  for (size_t i=0; i<matrix.rows()-1; i++) {
    Vector<Scalar> vec = matrix.row(i);
    stream << vec << ", ";
  }
  Vector<Scalar> vec = matrix.row(matrix.rows()-1);
  stream << vec << "]";

  return stream;
}

}
#endif // __LINALG_MATRIX_OPERATORS_HH__
