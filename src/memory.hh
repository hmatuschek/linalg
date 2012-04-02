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
 * A simple base class holding a owned or unowned reference to some data.
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
  virtual DataPtr<Scalar> *ref() = 0;

  virtual void unref() = 0;

  inline Scalar *ptr() {
    return _data;
  }
};


template<class Scalar>
class SmartPtr : public DataPtr
{
protected:
  bool _owns_data;
  size_t _refcount;

public:
  SmartPtr(Scalar *data, bool take_data)
    : _data(data), _owns_data(take_data), _refcount(1)
  {
    // Pass...
  }

  virtual DataPtr<Scalar> *ref()
  {
    this->_refcount++;
    return this;
  }

  virtual unref()
  {
    if (1 == _refcount) {
      _refcount=0;
      if (_owns_data)
        delete _data;
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
