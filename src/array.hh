#ifndef __LINALG_ARRAY_HH__
#define __LINALG_ARRAY_HH__

#include "memory.hh"
#include "exception.hh"
#include <vector>


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
protected:
  /** Manages the data. */
  DataPtr<Scalar> *_data;

  /** Offset of first element form _data ptr. */
  size_t _offset;

  /** Holds the dimensions of the array. */
  std::vector<size_t> _shape;

  /** Holds the strides of the array. */
  std::vector<size_t> _strides;


protected:
  /**
   * Hidden, full constructor.
   */
  Array(Scalar *data, size_t offset, const std::vector<size_t> &shape,
        const std::vector<size_t> &strides, bool take_data=true)
    : _data(new SmartPtr<Scalar>(data, take_data)), _offset(offset), _shape(shape), _strides(strides)
  {
    LINALG_SHAPE_ASSERT(shape.size() == strides.size());
  }


public:
  /**
   * Constructs uninitialized/empty array.
   */
  Array()
    : _data(0), _offset(0), _shape(0), _strides(0)
  {
    // Done..
  }


  /** Copy constructor. */
  Array(Array<Scalar> &other)
    : _data(0), _offset(other._offset), _shape(other._shape), _strides(other._strides)
  {
    if (0 != other._data) {
      _data = other._data->ref();
    }
  }


  /** Destructor. */
  virtual ~Array()
  {
    if (0 != _data) {
      _data->unref();
    }
  }


  /** Assignment operator. */
  Array<Scalar> &operator= (const Array<Scalar> &other)
  {
    // Unref "old" data, if there was some...
    if (0 != _data) {
      _data->unref();
    }

    // Copy ref and shape of other data...
    _data    = other._data;
    _shape   = other._shape;
    _strides = other._strides;

    return &this;
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
  inline Array<Scalar> t()
  {
    return Array<Scalar>(&this, _offset, std::vector<size_t>(_shape.rbegin(), _shape.rend()),
                         std::vector<size_t>(_strides.rbegin(), _strides.rend()));
  }
};


}

#endif // __LINALG_ARRAY_HH__
