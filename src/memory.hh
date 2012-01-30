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

  /**
   * If true, this instance owns the data and will free it if the instance gets destroyed.
   */
  bool   _owns_data;


public:
  /**
   * Constructs an empty data pointer (null pointer).
   */
  DataPtr() throw()
    : _data(0), _owns_data(false)
  {
    // Pass...
  }


  /**
   * constructor from data.
   */
  DataPtr(Scalar *data, bool owns_data) throw()
    : _data(data), _owns_data(owns_data)
  {
    // Pass...
  }


  /**
   * Constructs a weak reference to the data...
   */
  DataPtr(const DataPtr<Scalar> &ptr) throw()
    : _data(ptr._data), _owns_data(false)
  {
    // Pass...
  }


  /**
   * Copy constructor, that tries to take the ownership if @c taks_ownership is true.
   */
  DataPtr(DataPtr<Scalar> &ptr, bool take_ownership) throw(MemoryError)
    : _data(ptr._data), _owns_data(false)
  {
    if (take_ownership && ptr.ownsData()) {
      this->_owns_data = true;
      ptr.releaseData();
    }
    else if (take_ownership && ! ptr.ownsData()) {
      MemoryError err;
      err << "Can not take ownership of data: Source does not own the data.";
      throw err;
    }
  }


  /**
   * Allocates some new data.
   */
  explicit DataPtr(size_t size) throw()
    : _data(0), _owns_data(true)
  {
    this->_data = new Scalar[size];
  }


  /**
   * Frees the data if the instance holds the ownership of the data.
   */
  ~DataPtr() throw()
  {
    if (this->_owns_data && 0 != this->_data)
      delete this->_data;
  }


  /**
   * Returns true, if the instance holds the ownership of the data.
   */
  bool ownsData() const
  {
    return this->_owns_data;
  }


  /**
   * Releases the ownership of the data, does not free it!
   */
  void releaseData()
  {
    this->_data = 0;
    this->_owns_data = false;
  }


  /**
   * Returns true if there is no data.
   */
  bool isEmpty() const
  {
    return 0 == this->_data;
  }
};


}

#endif // __LINALG_MEMORY_HH__
