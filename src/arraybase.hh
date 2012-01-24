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
  size_t size;

  /**
   * Holds the pointer to the data owned by the container.
   */
  T* data;

  /**
   * If true, the array base owns the data and the data will freed if the instance gets destroyed.
   */
  bool owns_data;


public:
  ArrayBase(T *data, size_t size, bool owns_data)
    : data(data), size(size), owns_data(owns_data)
  {
    // Pass...
  }

  /**
   * Constructs a array with given size.
   */
  ArrayBase(size_t size)
    : size(size), data(0), owns_data(true)
  {
    // Just allocate the data
    this->data = new T[size];
  }

  /**
   * Copy constructor, takes the ownership, if the other array holds it.
   */
  ArrayBase(ArrayBase<T> &other)
    : size(other.size), data(other.data), owns_data(other.owns_data)
  {
    if (this->owns_data) {
      other.data = 0;
      other.owns_data = false;
    }
  }

  ArrayBase(const ArrayBase<T> &other)
    : size(other.size), data(other.data), owns_data(false)
  {
    // Pass...
  }

  /**
   * Destructor, also frees the data held.
   */
  ~ArrayBase()
  {
    if (this->owns_data)
      delete this->data;
  }

  /**
   * Assignment operator.
   */
  const ArrayBase<T> &operator= (ArrayBase<T> &other)
  {
    if (this->owns_data && 0 != this->data)
    {
      delete this->data;
    }

    this->owns_data = other.owns_data;
    this->data      = other.data;
    this->size      = other.size;

    if (this->owns_data) {
      other.owns_data = false;
      other.data = 0;
    }
  }

  /**
   * Returns the size of the array.
   */
  inline size_t getSize() const
  {
    return this->size;
  }

  /**
   * Returns a pointer to the data array.
   */
  inline T* operator* ()
  {
    return this->data;
  }

  /**
   * Retunrs a const pointer to the data array.
   */
  inline const T* operator* () const
  {
    return this->data;
  }

  /**
   * Retunrs the given element.
   */
  inline T &operator [](size_t i)
  {
    return this->data[i];
  }

  /**
   * Returns a weak reference to this array.
   */
  inline ArrayBase<T> weak() const
  {
    return ArrayBase<T>(this->data, this->size, false);
  }
};


}

#endif // ARRAYBASE_HH
