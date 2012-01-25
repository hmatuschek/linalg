/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

/**
 * @module mem Memory mangement
 *
 *
 */

#ifndef __LINALG_MEMORY_HH__
#define __LINALG_MEMORY_HH__

#include "exception.hh"
#include <memory>


namespace Linalg {

template <class Scalar>
class DataPtr
{
protected:
  Scalar *_data;
  bool   _owns_data;


public:
  DataPtr(Scalar *data, bool owns_data) throw()
    : _data(data), _owns_data(owns_data)
  {
    // Pass...
  }


  DataPtr(const DataPtr<Scalar> &ptr) throw()
    : _data(ptr._data), _owns_data(false)
  {
    // Pass...
  }


  DataPtr(DataPtr<Scalar> &ptr, bool take_ownership) throw()
    : _data(ptr._data), _owns_data(false)
  {
    if (take_ownership && ptr.ownsData())
    {
      this->_owns_data = true;
      ptr.releaseData();
    }
  }


  explicit DataPtr(size_t size) throw()
    : _data(0), _owns_data(true)
  {
    this->_data = new Scalar[size];
  }


  ~DataPtr() throw()
  {
    if (this->_owns_data && 0 != this->_data)
      delete this->_data;
  }


  bool ownsData() const
  {
    return this->_owns_data;
  }


  void releaseData()
  {
    this->_data = 0;
    this->_owns_data = false;
  }


  bool isEmpty() const
  {
    return 0 == this->_data;
  }
};


template <class MatrixType>
MatrixType takeOwnership(MatrixType &matrix)
{
  return MatrixType::unowned(new MatrixType(matrix, true));
}

}

#endif // __LINALG_MEMORY_HH__
