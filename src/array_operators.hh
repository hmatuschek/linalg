#ifndef __LINALG_ARRAY_OPERATORS_HH__
#define __LINALG_ARRAY_OPERATORS_HH__

#include "array.hh"

namespace Linalg {

/**
 * Comparison of arrays.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Array<bool> operator== (const Array<Scalar> &lhs, const Array<Scalar> &rhs)
{
  // Check dims:
  LINALG_SHAPE_ASSERT(lhs.ndims() == rhs.ndims());

  // Check shape:
  for (size_t i=0; i<lhs.ndims(); i++) {
    LINALG_SHAPE_ASSERT(lhs.dim(i) == rhs.dim(i));
  }

  // Allocate result array:
  Array<bool> res = res.empty(lhs.shape());

  // perform operation, element-wise:
  Array<Scalar>::iterator lhs_iter = lhs.begin();
  Array<Scalar>::iterator rhs_iter = rhs.begin();
  Array<bool>::iterator res_iter = res.begin();

  for (;! lhs_iter.atEnd(); lhs_iter++, rhs_iter++, res_iter++) {
    *res_iter = *lhs_iter == *rhs_iter;
  }

  return res;
}


bool any(const Array<bool> &array)
{
  Array<bool>::iterator iter = array.begin();
  for (; iter != array.end(); iter++) {
    if (*iter)
      return true;
  }
}


bool all(const Array<bool> &array)
{
  Array<bool>::iterator iter = array.begin();
  for (; iter != array.end(); iter++) {
    if (! *iter)
      return false;
  }
}


/**
 * General implementation of array sum.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Array<Scalar> operator+ (const Array<Scalar> &lhs, const Array<Scalar> &rhs)
{
  // Check dims:
  LINALG_SHAPE_ASSERT(lhs.ndims() == rhs.ndims());

  // Check shape:
  for (size_t i=0; i<lhs.ndims(); i++) {
    LINALG_SHAPE_ASSERT(lhs.dim(i) == rhs.dim(i));
  }

  // Allocate result array:
  Array<Scalar> res = res.empty(lhs.shape());

  // perform operation, element-wise:
  Array<Scalar>::iterator lhs_iter = lhs.begin();
  Array<Scalar>::iterator rhs_iter = rhs.begin();
  Array<Scalar>::iterator res_iter = res.begin();

  for (;! lhs_iter.atEnd(); lhs_iter++, rhs_iter++, res_iter++) {
    *res_iter = *lhs_iter + *rhs_iter;
  }

  return res;
}


/**
 * General implementation of array inplace-sum.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Array<Scalar> &operator+= (Array<Scalar> &lhs, const Array<Scalar> &rhs)
{
  // Check dims:
  LINALG_SHAPE_ASSERT(lhs.ndims() == rhs.ndims());

  // Check shape:
  for (size_t i=0; i<lhs.ndims(); i++) {
    LINALG_SHAPE_ASSERT(lhs.dim(i) == rhs.dim(i));
  }

  // perform operation, element-wise:
  Array<Scalar>::iterator lhs_iter = lhs.begin();
  Array<Scalar>::iterator rhs_iter = rhs.begin();

  for (;! lhs_iter.atEnd(); lhs_iter++, rhs_iter++) {
    *lhs_iter += *rhs_iter;
  }

  return res;
}

}


#endif // __LINALG_ARRAY_OPERATORS_HH__
