#ifndef __LINALG_ARRAY_HH__
#define __LINALG_ARRAY_HH__

#include "memory.hh"
#include "exception.hh"
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
class Array
{
public:
  /**
   * Iterator class, to iterate over all elements of an array.
   */
  class iterator {
  private:
    /** Holds the current index. */
    std::vector<size_t> _current_idx;

    /** Holds a reference to the array. */
    Array<Scalar> &_array;

  public:
    /**
     * Constructor.
     */
    iterator(Array<Scalar> &array, const std::vector<size_t> &idx)
      : _current_idx(idx), _array(array)
    {
      // Pass...
    }

    /**
     * Assignment operator.
     */
    iterator &operator= (iterator &other)
    {
      _current_idx     = other._current_idx;
      _array           = other._array;

      return *this;
    }

    /**
     * Increment operator, to shif the iterator to the next element.
     */
    iterator &operator++ (int)
    {
      // Check if iterator is at the end:
      if(_current_idx.back() == _array.shape().back())
        return *this;

      // Otherwise:
      _current_idx[0]++;
      for(size_t i=0; i<_array.ndim()-1; i++) {
        if (_current_idx[i]==_array.shape(i)) {
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
    inline bool operator==(const iterator &other)
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
    bool operator!= (const iterator &other)
    {
      return ! (*this == other);
    }

    /**
     * Dereferencing the element addressed by the iterator.
     */
    Scalar &operator*()
    {
      size_t offset = 0;
      for (size_t i=0; i<_array.ndim(); i++) {
        offset += _current_idx[i]*_array.strides(i);
      }

      return *(_array.ptr() + offset);
    }
  };

  /**
   * Const iterator class, to iterate over all elements of an array.
   */
  class const_iterator {
  private:
    /** Holds the current index. */
    std::vector<size_t> _current_idx;

    /** Holds a reference to the array. */
    const Array<Scalar> &_array;

  public:
    /**
     * Constructor.
     */
    const_iterator(const Array<Scalar> &array, const std::vector<size_t> &idx)
      : _current_idx(idx), _array(array)
    {
      // Pass...
    }

    /**
     * Assignment operator.
     */
    const_iterator &operator= (const const_iterator &other)
    {
      _current_idx     = other._current_idx;
      _array           = other._array;

      return *this;
    }

    /**
     * Increment operator, to shif the iterator to the next element.
     */
    const_iterator &operator++ (int)
    {
      // Check if iterator is at the end:
      if(_current_idx.back() == _array.shape().back())
        return *this;

      // Otherwise:
      _current_idx[0]++;
      for(size_t i=0; i<_array.ndim()-1; i++) {
        if (_current_idx[i]==_array.shape(i)) {
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
    inline bool operator==(const const_iterator &other)
    {
      return _current_idx == other._current_idx
          && _array.data()->ptr() == other._array.data()->ptr()
          && _array.offset() == other._array.offset()
          && _array.shape() == other._array.shape()
          && _array.strides() == other._array.strides();
    }

    /**
     * Iterator comparison.
     */
    bool operator!= (const const_iterator &other)
    {
      return ! (*this == other);
    }

    /**
     * Dereferencing the element addressed by the iterator.
     */
    const Scalar &operator*() const
    {
      size_t offset = 0;
      for (size_t i=0; i<_array.ndim(); i++) {
        offset += _current_idx[i]*_array.strides(i);
      }
      return *(_array.ptr() + offset);
    }
  };

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
  /** Manages the data. */
  DataPtr<Scalar> *_data;

  /** Offset of first element form _data ptr. */
  size_t _offset;

  /** Holds the dimensions of the array. */
  std::vector<size_t> _shape;

  /** Holds the strides of the array. */
  std::vector<size_t> _strides;

  /**
   * Holds and value-reference instance to the data.
   */
  Values _values;


protected:
  /**
   * Hidden, full constructor from pointer.
   */
  Array(Scalar *data, size_t offset, const std::vector<size_t> &shape,
        const std::vector<size_t> &strides, bool take_data)
    : _data(new SmartPtr<Scalar>(data, take_data)), _offset(offset),
      _shape(shape), _strides(strides), _values(*this)
  {
    LINALG_SHAPE_ASSERT(shape.size() == strides.size());
  }


public:
  /**
   * Full constructor from data.
   */
  Array(DataPtr<Scalar> *data, size_t offset, const std::vector<size_t> &shape,
        const std::vector<size_t> &strides)
    : _data(data->ref()), _offset(offset), _shape(shape), _strides(strides),
      _values(*this)
  {
    LINALG_SHAPE_ASSERT(shape.size() == strides.size());
  }


  /**
   * Constructs uninitialized/empty array.
   */
  Array()
    : _data(0), _offset(0), _shape(0), _strides(0), _values(*this)
  {
    // Done..
  }


  Array(const std::vector<size_t> &shape, bool rowmajor=true)
    : _data(0), _offset(0), _shape(shape), _strides(shape.size()), _values(*this)
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

    // Allocate some data:
    _data = new SmartPtr<Scalar>(new Scalar[size], true);

    // Done...
  }

  /** Copy constructor. */
  Array(const Array<Scalar> &other)
    : _data(other._data->ref()), _offset(other._offset), _shape(other._shape), _strides(other._strides),
      _values(*this)
  {
    // Pass...
  }


  /** Destructor. */
  virtual ~Array()
  {
    if (0 != _data) {
      _data->unref();
      this->_data=0;
    }
  }


  /** Assignment operator. */
  Array<Scalar> &operator= (const Array<Scalar> &other)
  {
    // Get reference to the data of other:
    DataPtr<Scalar> *_new_data = other._data->ref();

    // Unref "old" data, if there was some...
    if (0 != _data) {
      _data->unref();
    }

    // Copy ref and shape of other data...
    this->_data    = _new_data;
    this->_shape   = other._shape;
    this->_strides = other._strides;
    this->_offset  = other._offset;

    return *this;
  }


  Array<Scalar> copy()
  {
    Array cp(_shape);
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


  /** Returns the element at the given index. */
  Scalar &at(const std::vector<size_t> &idxs)
  {
    LINALG_SHAPE_ASSERT(idxs.size() == _strides.size());

    size_t offset = _offset;
    for (size_t i=0; i<_shape.size(); i++) {
      offset += idxs[i]*_strides[i];
    }

    return _data->ptr()[offset];
  }


  /** Returns the element at the given index. */
  const Scalar &at(const std::vector<size_t> &idxs) const
  {
    LINALG_SHAPE_ASSERT(idxs.size() == _strides.size());

    size_t offset = _offset;
    for (size_t i=0; i<_shape.size(); i++) {
      offset += idxs[i]*_strides[i];
    }

    return _data->prt()[offset];
  }


  /** Returns the transpose of the array. */
  inline Array<Scalar> t() const
  {
    return Array<Scalar>(this->_data, this->_offset,
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
   */
  inline Values &operator* () {
    return _values;
  }


  /**
   * Returns the pointer to the first element.
   */
  inline const Scalar *ptr() const {
    return this->_data->ptr() + this->_offset;
  }

  /**
   * Returns the pointer to the first element.
   */
  inline Scalar *ptr() {
    return this->_data->ptr() + this->_offset;
  }

  inline iterator begin() {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    return iterator(*this, idxs);
  }

  inline iterator end() {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    idxs.back() = this->shape().back();
    return iterator(*this, idxs);
  }

  inline const_iterator const_begin() const {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    return const_iterator(*this, idxs);
  }

  inline const_iterator const_end() const {
    std::vector<size_t> idxs(this->_shape.size(), 0);
    idxs.back() = this->shape().back();
    return const_iterator(*this, idxs);
  }
};

}

#endif // __LINALG_ARRAY_HH__
