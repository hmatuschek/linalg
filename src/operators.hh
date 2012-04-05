/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */


#ifndef __LINALG_OPERATORS_HH__
#define __LINALG_OPERATORS_HH__

#include "matrix.hh"
#include "blas/gemm.hh"
#include "array_operators.hh"
#include "trimatrix_operators.hh"


namespace Linalg {


/**
 * Implements the common matrix product.
 *
 * @ingroup operators
 */
inline Matrix<double>
operator* (const Matrix<double> &lhs, const Matrix<double> &rhs)
throw (ShapeError)
{
  // Allocate matrix for result:
  Matrix<double> result(lhs.rows(), rhs.cols());
  // Use Blas::gemm() to compute product
  Blas::gemm(1., lhs, rhs, 0.0, result);
  // Pass ownership of result-matrix to caller...
  return result;
}


/**
 * Implements the common matrix sum.
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


template <class Scalar>
std::ostream &operator<< (std::ostream &stream, Vector<Scalar> &vector)
{
  stream << "[";
  for (size_t i=0; i<vector.dim()-1; i++) {
    stream << vector(i) << ", ";
  }
  stream << vector(vector.dim()-1) << "]";

  return stream;
}


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

#endif // __LINALG_OPERATORS_HH__
