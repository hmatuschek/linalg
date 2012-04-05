#ifndef __LINALG_TRIMATRIX_OPERATORS_HH__
#define __LINALG_TRIMATRIX_OPERATORS_HH__

#include "trimatrix.hh"

namespace Linalg {

/**
 * Simple, element wise comparison.
 *
 * @bug Does not work properly.
 */
template <class Scalar>
Matrix<bool> operator== (const TriMatrix<Scalar> &lhs, const TriMatrix<Scalar> &rhs)
{
  LINALG_SHAPE_ASSERT(lhs.ndim() == rhs.ndim());
  for (size_t i=0; i<lhs.ndim(); i++) {
    LINALG_SHAPE_ASSERT(lhs.shape(i) == rhs.shape(i));
  }

  Matrix<bool> res(lhs.rows(), lhs.cols());

  for (size_t i=0; i<lhs.rows(); i++) {
    for (size_t j=0; j<lhs.cols(); j++) {
      res(i,j) = (lhs(i,j) == rhs(i,j));
    }
  }

  return res;
}


}


#endif // __LINALG_TRIMATRIX_OPERATORS_HH__
