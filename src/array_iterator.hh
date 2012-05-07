/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_ARRAY_ITERATOR_HH__
#define __LINALG_ARRAY_ITERATOR_HH__

#include <vector>


namespace Linalg {
/* Forward declaration. */
template<class Scalar> class Array;


/**
 * Implements iteration over complete array.
 */
template <class Scalar>
class ArrayIterator
{
private:
  /** Holds the current index. */
  std::vector<size_t> _current_idx;

  /** Holds a reference to the array. */
  Array<Scalar> *_array;

public:
  /** Default constructor. */
  ArrayIterator()
    : _current_idx(0), _array(0)
  {
    // Pass...
  }

  /**
   * Constructor.
   */
  ArrayIterator(Array<Scalar> *array, const std::vector<size_t> &idx)
    : _current_idx(idx), _array(array)
  {
    // Pass...
  }

  /**
   * Assignment operator.
   */
  ArrayIterator<Scalar> &operator= (ArrayIterator<Scalar> &other)
  {
    _current_idx     = other._current_idx;
    _array           = other._array;

    return *this;
  }

  /**
   * Increment operator, to shif the iterator to the next element.
   */
  ArrayIterator<Scalar> &operator++ (int)
  {
    // Check if iterator is at the end:
    if(_current_idx.back() == _array->shape().back())
      return *this;

    // Otherwise:
    _current_idx[0]++;
    for(size_t i=0; i<_array->ndim()-1; i++) {
      if (_current_idx[i]==_array->shape(i)) {
        _current_idx[i] = 0;
        _current_idx[i+1]++;
      }
    }

    // Done.
    return *this;
  }

  /**
   * Iterator comparison.
   */
  inline bool operator==(const ArrayIterator<Scalar> &other)
  {
    if (_current_idx.size() != other._current_idx.size())
      return false;
    for (size_t i=0; i<_current_idx.size(); i++) {
      if (_current_idx[i] != other._current_idx[i])
        return false;
    }

    return true;
  }

  /**
   * Iterator comparison.
   */
  bool operator!= (const ArrayIterator<Scalar> &other)
  {
    return ! (*this == other);
  }

  /**
   * Dereferencing the element addressed by the iterator.
   */
  Scalar &operator*()
  {
    size_t offset = 0;
    for (size_t i=0; i<_array->ndim(); i++) {
      offset += _current_idx[i]*_array->strides(i);
    }

    return *(_array->ptr() + offset);
  }
};



/**
 * Const iterator class, to iterate over all elements of an array.
 */
template <class Scalar>
class ArrayConstIterator {
private:
  /** Holds the current index. */
  std::vector<size_t> _current_idx;

  /** Holds a reference to the array. */
  const Array<Scalar> *_array;

public:
  /** Default constructor. */
  ArrayConstIterator() : _current_idx(0), _array(0) { /* Pass... */ }

  /**
   * Constructor.
   */
  ArrayConstIterator(const Array<Scalar> *array, const std::vector<size_t> &idx)
    : _current_idx(idx), _array(array)
  {
    // Pass...
  }

  /**
   * Assignment operator.
   */
  ArrayConstIterator &operator= (const ArrayConstIterator &other)
  {
    _current_idx     = other._current_idx;
    _array           = other._array;

    return *this;
  }

  /**
   * Increment operator, to shif the iterator to the next element.
   */
  ArrayConstIterator &operator++ (int)
  {
    // Check if iterator is at the end:
    if(_current_idx.back() == _array->shape().back())
      return *this;

    // Otherwise:
    _current_idx[0]++;
    for(size_t i=0; i<_array->ndim()-1; i++) {
      if (_current_idx[i]==_array->shape(i)) {
        _current_idx[i] = 0;
        _current_idx[i+1]++;
      }
    }

    // Done.
    return *this;
  }

  /**
   * Iterator comparison.
   */
  inline bool operator==(const ArrayConstIterator &other)
  {
    if (_current_idx.size() != other._current_idx.size())
      return false;
    for (size_t i=0; i<_current_idx.size(); i++) {
      if (_current_idx[i] != other._current_idx[i])
        return false;
    }

    return true;
  }

  /**
   * Iterator comparison.
   */
  bool operator!= (const ArrayConstIterator &other)
  {
    return ! (*this == other);
  }

  /**
   * Dereferencing the element addressed by the iterator.
   */
  const Scalar &operator*() const
  {
    size_t offset = 0;
    for (size_t i=0; i<_array->ndim(); i++) {
      offset += _current_idx[i]*_array->strides(i);
    }
    return *(_array->ptr() + offset);
  }
};

}


#endif // ARRAY_ITERATOR_HH
