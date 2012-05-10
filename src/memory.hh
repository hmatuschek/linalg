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
 * Template prototype for manager classes.
 */
template <class Scalar>
class DataMngr
{
private:
  Scalar *_data;
  bool _owned;

public:
  DataMngr(Scalar *data, bool owned)
    : _data(data), _owned(owned)
  {
    // Pass...
  }

  explicit DataMngr(size_t size)
  {
    _data = new Scalar[size];
    _owned = true;
  }

  virtual ~DataMngr() {
    if (_owned) {
      delete _data;
    }
  }

  Scalar *ptr() const {
    return _data;
  }

  bool ownsData() const {
    return _owned;
  }
};



/**
 * A simple base class holding an owned or unowned reference to some data.
 */
template <class Scalar>
class DataPtr
{
protected:
  /**
   * Holds a reference to the data.
   */
  DataMngr<Scalar> *_mngr;

  /**
   * Holds the actual data pointer.
   */
  Scalar *_data;

  /**
   * Holds a reference to the reference counter.
   */
  size_t *_refcount;


public:
  /**
   * Empty constructor.
   */
  DataPtr()
    : _mngr(0), _data(0), _refcount(0)
  {
    // Pass...
  }

  /**
   * Constructor.
   */
  DataPtr(DataMngr<Scalar> *mngr)
    : _mngr(mngr), _data(mngr->ptr())
  {
    _refcount = new size_t(1);
  }

  /**
   * Copy constructor.
   */
  DataPtr(const DataPtr<Scalar> &other)
    : _mngr(other._mngr), _data(other._data), _refcount(other._refcount)
  {
    if (0 != _refcount)
      (*_refcount)++;
  }

  /**
   * Destructor.
   */
  ~DataPtr()
  {
    if (0 != _refcount) {
      (*_refcount)--;
      if (0 == *_refcount) {
        delete _mngr;
        delete _refcount;
        _mngr = 0;
        _refcount = 0;
        _data = 0;
      }
    }
  }


  /**
   * Assignment operator.
   */
  DataPtr<Scalar> &operator =(const DataPtr<Scalar> &other)
  {
    if (0 != _refcount) {
      (*_refcount)--;
      if (0 == *_refcount) {
        delete _mngr;
        delete _refcount;
        _mngr = 0;
        _refcount = 0;
        _data = 0;
      }
    }

    _mngr = other._mngr;
    _data = other._data;
    _refcount = other._refcount;
    if (0 != _refcount)
      (*_refcount)++;

    return *this;
  }
};

}

#endif // __LINALG_MEMORY_HH__
