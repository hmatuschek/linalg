/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_ARRAY_HH__
#define __LINALG_ARRAY_HH__

#include "memory.hh"
#include "exception.hh"
#include "array_iterator.hh"
#include <vector>

#include <iostream>


namespace Linalg {


/**
 * This class represents an N-dimensional array.
 *
 * The implementation is very close to the ndarray as implemented in NumPy, this allows for an
 * easy encapsulation of NumPy array in this template library.
 */
template <class Scalar>
class Array : public DataPtr<Scalar>
{
public:
  /**
   * Iterator class, to iterate over all elements of an array.
   */
  typedef ArrayIterator<Scalar> iterator;

  /**
   * Const iterator class, to iterate over all elements of an array.
   */
  typedef ArrayConstIterator<Scalar> const_iterator;

  /**
   * This class represents the values of the array.
   */
  class Values {
  protected:
    /** Holds a weak reference to the array. */
    Array<Scalar> &_array;

  public:
    /** Constructor. */
    Values(Array<Scalar> &array): _array(array) { /* Pass... */ }

    /** Assignment. */
    Values &operator= (const Array<Scalar> &other)
    {
      // Check array dimension:
      LINALG_SHAPE_ASSERT(_array.ndim() == other.ndim());
      for (size_t i=0; i<_array.ndim(); i++) {
        LINALG_SHAPE_ASSERT(_array.shape(i) == other.shape(i));
      }

      // Iterate over all values of _array and array and assign other -> _array:
      iterator a_iter = _array.begin();
      const_iterator b_iter = other.const_begin();
      for (; a_iter != _array.end(); a_iter++, b_iter++) {
        *a_iter = *b_iter;
      }

      // Done..
      return *this;
    }

    /** Assignment. */
    Values &operator= (const Scalar &value) {
      // Iterate over all values of _array assign value -> _array:
      iterator a_iter = _array.begin();
      for (; a_iter != _array.end(); a_iter++) {
        *a_iter = value;
      }

      // Done..
      return *this;
    }
  };


protected:
  /** Offset of first element form _data ptr. */
  size_t _offset;

  /** Holds the dimensions of the array. */
  std::vector<size_t> _shape;

  /** Holds the strides of the array. */
  std::vector<size_t> _strides;

  /**
   * Holds a value-reference instance to the data.
   */
  Values _values;


protected:
  /**
   * Hidden, full constructor from pointer.
   */
  Array(Scalar *data, size_t offset, const std::vector<size_t> &shape,
        const std::vector<size_t> &strides, bool take_data)
    : DataPtr<Scalar>(new DataMngr<Scalar>(data, take_data)), _offset(offset),
      _shape(shape), _strides(strides), _values(*this)
  {
    LINALG_SHAPE_ASSERT(shape.size() == strides.size());
  }


public:
  /**
   * Full constructor from data.
   */
  Array(const DataPtr<Scalar> &data, size_t offset, const std::vector<size_t> &shape,
        const std::vector<size_t> &strides)
    : DataPtr<Scalar>(data), _offset(offset), _shape(shape), _strides(strides),
      _values(*this)
  {
    LINALG_SHAPE_ASSERT(shape.size() == strides.size());
  }


  /**
   * Constructs uninitialized/empty array.
   */
  Array()
    : DataPtr<Scalar>(), _offset(0), _shape(0), _strides(0), _values(*this)
  {
    // Done..
  }


  /**
   * Allocates a new (empty) array of the given shape.
   */
  Array(const std::vector<size_t> &shape, bool rowmajor=true)
    : DataPtr<Scalar>(), _offset(0), _shape(shape), _strides(shape.size()), _values(*this)
  {
    // Check if dimension
    if (0 == _shape.size())
      return;

    // Determine size and create strides in column-major:
    size_t size = 1;
    if (rowmajor) {
      std::vector<size_t>::reverse_iterator sh_iter = _shape.rbegin();
      std::vector<size_t>::reverse_iterator st_iter = _strides.rbegin();
      for(; sh_iter != _shape.rend(); sh_iter++, st_iter++) {
        *st_iter = size;
        size *= *sh_iter;
      }
    } else {
      std::vector<size_t>::iterator sh_iter = _shape.begin();
      std::vector<size_t>::iterator st_iter = _strides.begin();
      for(; sh_iter != _shape.end(); sh_iter++, st_iter++) {
        *st_iter = size;
        size *= *sh_iter;
      }
    }

    // Allocate some data and assign it to this:
    DataPtr<Scalar>::operator =(DataPtr<Scalar>(new DataMngr<Scalar>(size)));

    // Done...
  }


  /** Copy constructor. */
  Array(const Array<Scalar> &other)
    : DataPtr<Scalar>(other), _offset(other._offset), _shape(other._shape), _strides(other._strides),
      _values(*this)
  {
    // Pass...
  }


  /**
   * Assignment operator.
   *
   * This operator does not copy any values from the given array into this array. Arrays are
   * containers holding the elements and some meta-data about the shape and storage of the elements.
   *
   * Therefore, after assignment, both array share the same elements and have identical shape. This
   * is a array view assignment.
   *
   * To copy values, use either @c copy to create a new copy of the array, or @c values to
   * set all values of this array at once.
   */
  Array<Scalar> &operator= (const Array<Scalar> &other)
  {
    // Get reference to the data of other:
    DataPtr<Scalar>::operator =(other);

    // Copy ref and shape of other data...
    this->_shape   = other._shape;
    this->_strides = other._strides;
    this->_offset  = other._offset;

    return *this;
  }


  /**
   * Creates a copy of this array.
   */
  Array<Scalar> copy(bool rowmajor=true)
  {
    Array cp(_shape, rowmajor);
    cp.values() = *this;
    return cp;
  }


  /**
   * Retruns the number of dimensions of the array.
   */
  inline size_t ndim() const
  {
    return _shape.size();
  }


  /**
   * Returns the shape of the array.
   */
  inline const std::vector<size_t> &shape() const
  {
    return _shape;
  }


  /**
   * Returns the size of the i-th dimension.
   */
  inline size_t shape(size_t i) const
  {
    return _shape[i];
  }


  /**
   * Returns the offset of the first element from the internal pointer.
   */
  inline size_t offset() const
  {
    return _offset;
  }


  /**
   * Returns the strides of the array.
   */
  inline const std::vector<size_t> &strides() const
  {
    return _strides;
  }


  /**
   * Returns the i-th stride.
   */
  inline size_t strides(size_t i) const
  {
    return _strides[i];
  }


  /**
   * Returns true, if the storage-order of the array is row-major (C order).
   */
  inline bool isRowMajor() const {
    return 1 == this->_strides.back();
  }


  /**
   * Retunrs the managed data pointer.
   */
  inline DataPtr<Scalar> *data() const
  {
    return this->_data;
  }


  /**
   * Returns the element at the given index.
   */
  Scalar &at(const std::vector<size_t> &idxs)
  {
    LINALG_SHAPE_ASSERT(idxs.size() == _strides.size());

    size_t offset = _offset;
    for (size_t i=0; i<_shape.size(); i++) {
      offset += idxs[i]*_strides[i];
    }

    return this->_data[offset];
  }


  /**
   * Returns the element at the given index-vector.
   */
  const Scalar &at(const std::vector<size_t> &idxs) const
  {
    LINALG_SHAPE_ASSERT(idxs.size() == _strides.size());

    size_t offset = _offset;
    for (size_t i=0; i<_shape.size(); i++) {
      offset += idxs[i]*_strides[i];
    }

    return this->_data->prt()[offset];
  }


  /**
   * Returns the transposed of the array.
   */
  inline Array<Scalar> t() const
  {
    return Array<Scalar>(*this, this->_offset,
                         std::vector<size_t>(_shape.rbegin(), _shape.rend()),
                         std::vector<size_t>(_strides.rbegin(), _strides.rend()));
  }


  /**
   * Returns a reference to the values of the array.
   */
  inline Values &values() {
    return _values;
  }

  /**
   * Returns a reference to the values of the array.
   *
   * Equivalent to call @c values.
   */
  inline Values &operator* () {
    return _values;
  }


  /**
   * Returns the pointer to the first element.
   */
  inline const Scalar *ptr() const {
    return this->_data + this->_offset;
  }

  /**
   * Returns the pointer to the first element.
   */
  inline Scalar *ptr() {
    return this->_data + this->_offset;
  }


  /**
   * Returns an iterator over all elements of the array, pointing to the first element (0,0,...).
   */
  inline iterator begin() {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    return iterator(this, idxs);
  }

  /**
   * Returns an iterator over all elements of the array, pointing right after the last element
   * of the array.
   */
  inline iterator end() {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    idxs.back() = this->shape().back();
    return iterator(this, idxs);
  }

  /**
   * Returns a const iterator over all elements pointing to the first element (0,0,...).
   */
  inline const_iterator const_begin() const {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    return const_iterator(this, idxs);
  }

  /**
   * Returns a const iterator over all elements of the array pointing right after the last element.
   */
  inline const_iterator const_end() const {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    idxs.back() = this->shape().back();
    return const_iterator(this, idxs);
  }
};

}

#endif // __LINALG_ARRAY_HH__
