/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_ARRAY_OPERATORS_HH__
#define __LINALG_ARRAY_OPERATORS_HH__

#include "array.hh"


namespace Linalg {

/**
 * Comparison of arrays.
 *
 * @ingroup operators
 */
template <class T>
Array<bool> operator== (const Array<T> &lhs, const Array<T> &rhs)
{
  // Check dims:
  LINALG_SHAPE_ASSERT(lhs.ndim() == rhs.ndim());

  // Check shape:
  for (size_t i=0; i<lhs.ndim(); i++) {
    LINALG_SHAPE_ASSERT(lhs.shape(i) == rhs.shape(i));
  }

  // Allocate result array:
  Array<bool> res(lhs.shape());

  // perform operation, element-wise:
  ArrayConstIterator<T> lhs_iter = lhs.const_begin();
  ArrayConstIterator<T> rhs_iter = rhs.const_begin();
  Array<bool>::iterator res_iter = res.begin();

  for (; lhs_iter != lhs.const_end(); lhs_iter++, rhs_iter++, res_iter++) {
    *res_iter = (*lhs_iter == *rhs_iter);
  }

  return res;
}


/**
 * This function returns true if at least one elements of the given array of booleans is true.
 */
inline bool any(const Array<bool> &array)
{
  Array<bool>::const_iterator iter = array.const_begin();
  for (; iter != array.const_end(); iter++) {
    if (*iter)
      return true;
  }

  return false;
}


/**
 * This function retunrs true if all elements of the given array of booleans are true.
 */
inline bool all(const Array<bool> &array)
{
  Array<bool>::const_iterator iter = array.const_begin();
  for (; iter != array.const_end(); iter++) {
    if (! *iter)
      return false;
  }

  return true;
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
  Array<Scalar> res(lhs.shape());

  // perform operation, element-wise:
  ArrayConstIterator<Scalar> lhs_iter = lhs.const_begin();
  ArrayConstIterator<Scalar> rhs_iter = rhs.const_begin();
  ArrayIterator<Scalar> res_iter = res.begin();

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
template <class T>
inline Array<T> &operator+= (Array<T> &lhs, const Array<T> &rhs)
{
  // Check dims:
  LINALG_SHAPE_ASSERT(lhs.ndims() == rhs.ndims());

  // Check shape:
  for (size_t i=0; i<lhs.ndims(); i++) {
    LINALG_SHAPE_ASSERT(lhs.dim(i) == rhs.dim(i));
  }

  // perform operation, element-wise:
  ArrayIterator<T> lhs_iter = lhs.begin();
  ArrayConstIterator<T> rhs_iter = rhs.const_begin();

  for (;! lhs_iter.atEnd(); lhs_iter++, rhs_iter++) {
    *lhs_iter += *rhs_iter;
  }

  return lhs;
}

}


#endif // __LINALG_ARRAY_OPERATORS_HH__
