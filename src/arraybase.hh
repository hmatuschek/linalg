/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */


#ifndef __LINALG_ARRAYBASE_HH__
#define __LINALG_ARRAYBASE_HH__

#include <cstdlib>
#include <cstring>


namespace Linalg {


/**
 * Just a container holding the data as a dense 1D array.
 *
 * @ingroup core
 */
template<class T>
class ArrayBase
{
protected:
  /**
   * Holds the size of the data array.
   */
  size_t _size;

  /**
   * Holds the pointer to the data owned by the container.
   */
  T* _data;

  /**
   * If true, the array base owns the data and the data will freed if the instance gets destroyed.
   */
  bool _owns_data;


public:
  ArrayBase(T *data, size_t size, bool owns_data)
    : _data(data), _size(size), _owns_data(owns_data)
  {
    // Pass...
  }

  /**
   * Constructs a array with given size.
   */
  ArrayBase(size_t size)
    : _size(size), _data(0), _owns_data(true)
  {
    // Just allocate the data
    this->_data = new T[size];
  }

  /**
   * Copy constructor, takes the ownership, if the other array holds it.
   */
  ArrayBase(ArrayBase<T> &other)
    : _size(other._size), _data(other._data), _owns_data(other._owns_data)
  {
    if (this->_owns_data) {
      other._data = 0;
      other._owns_data = false;
    }
  }

  ArrayBase(const ArrayBase<T> &other)
    : _size(other._size), _data(other._data), _owns_data(false)
  {
    // Pass...
  }

  /**
   * Destructor, also frees the data held.
   */
  ~ArrayBase()
  {
    if (this->_owns_data)
      delete this->_data;
  }

  /**
   * Assignment operator.
   */
  const ArrayBase<T> &operator= (ArrayBase<T> &other)
  {
    if (this->_owns_data && 0 != this->_data)
    {
      delete this->_data;
    }

    this->_owns_data = other._owns_data;
    this->_data      = other._data;
    this->_size      = other._size;

    if (this->_owns_data) {
      other._owns_data = false;
      other._data = 0;
    }
  }


  /**
   * Assignment (weak) operator.
   */
  const ArrayBase<T> &operator= (const ArrayBase<T> &other)
  {
    if (this->_owns_data && 0 != this->_data)
    {
      delete this->_data;
    }

    this->_owns_data = false;
    this->_data      = other._data;
    this->_size      = other._size;
  }


  /**
   * Returns the size of the array.
   */
  inline size_t getSize() const
  {
    return this->_size;
  }

  /**
   * Returns a pointer to the data array.
   */
  inline T* operator* ()
  {
    return this->_data;
  }

  /**
   * Retunrs a const pointer to the data array.
   */
  inline const T* operator* () const
  {
    return this->_data;
  }

  /**
   * Retunrs the given element.
   */
  inline T &operator [](size_t i)
  {
    return this->_data[i];
  }

  /**
   * Returns a weak reference to this array.
   */
  inline ArrayBase<T> weak() const
  {
    return ArrayBase<T>(this->_data, this->_size, false);
  }
};


}

#endif // ARRAYBASE_HH
