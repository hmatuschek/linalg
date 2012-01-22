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


public:
  /**
   * Constructs a array with given size.
   */
  ArrayBase(size_t size)
    : size(size), data(0)
  {
    // Just allocate the data
    this->data = new T[size];
  }

  /**
   * Copy constructor, also copies the content of the given array.
   */
  ArrayBase(const ArrayBase<T> &other)
    : size(other.size), data(0)
  {
    this->data = new T[this->size];
    memcpy(this->data, other.data, sizeof(T)*this->size);
  }

  /**
   * Destructor, also frees the data held.
   */
  ~ArrayBase()
  {
    delete this->data;
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
  inline T *getPtr()
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
};


}

#endif // ARRAYBASE_HH
