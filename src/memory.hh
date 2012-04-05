/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */


#ifndef __LINALG_MEMORY_HH__
#define __LINALG_MEMORY_HH__

#include "exception.hh"
#include <memory>


namespace Linalg {

/**
 * A simple base class holding an owned or unowned reference to some data.
 */
template <class Scalar>
class DataPtr
{
protected:
  /**
   * Holds the actual data pointer.
   */
  Scalar *_data;


public:
  DataPtr(Scalar *ptr)
    : _data(ptr)
  {
    // Pass...
  }

  virtual ~DataPtr() { /* pass... */ }

  virtual DataPtr<Scalar> *ref() = 0;

  virtual void unref() = 0;

  inline Scalar *ptr() {
    return _data;
  }
};



/**
 * Simple reference counting data pointer.
 */
template<class Scalar>
class SmartPtr : public DataPtr<Scalar>
{
protected:
  /** If true, the data will be freed if the reference count gets 0. */
  bool _owns_data;
  /** The reference counter. */
  size_t _refcount;

public:
  /**
   * Constructor.
   *
   * @param data Specifies the managed data pointer.
   * @param take_data Specifies if the data should be freed if not used anymore.
   */
  SmartPtr(Scalar *data, bool take_data)
    : DataPtr<Scalar>(data), _owns_data(take_data), _refcount(1)
  {
    // Pass...
  }

  /**
   * Creates a "new" reference to the data.
   */
  virtual DataPtr<Scalar> *ref()
  {
    _refcount++;
    return this;
  }

  /**
   * Dereferences the data.
   */
  virtual void unref()
  {
    if (1 == _refcount) {
      _refcount=0;
      if (_owns_data)
        delete[] this->ptr();
      delete this;
    }
    else
    {
      _refcount--;
    }
  }
};

}

#endif // __LINALG_MEMORY_HH__
