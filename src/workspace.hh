/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_WORKSPACE_HH__
#define __LINALG_WORKSPACE_HH__

#include <cstdlib>


namespace Linalg {


/**
 * Semi-automatic workspace allocation for LAPACK routines.
 *
 * @ingroup linalg
 */
class Workspace
{
protected:
  /**
   * Holds the size of the allocated memory (in bytes).
   */
  size_t _size;

  /**
   * Holds a pointer to the allocated workspace memory.
   */
  char *_data;


public:
  /**
   * Creates a workspace with the given size (in bytes).
   */
  Workspace(size_t size=0)
    : _size(size), _data(0)
  {
    // Preallocate some space:
    if (0 < size)
    {
      this->_data = new char[size];
    }
  }

  /**
   * Frees the allocated memory.
   */
  ~Workspace()
  {
    // Free workspace if allocated:
    if (0 != this->_data)
    {
      delete[] this->_data;
    }
  }

  /**
   * Ensures that at least @c size elements are allocated for the working-memory.
   */
  template <class Scalar>
  Scalar *ensure(size_t size)
  {
    // If more space is requested than allocated -> allocate more
    if (this->_size < size*sizeof(Scalar))
    {
      if (0 != this->_data)
      {
        delete[] this->_data;
      }

      this->_data = (void *)new Scalar[size];
      this->_size = size*sizeof(Scalar);
    }

    return (Scalar *)_data;
  }

  template <class Scalar>
  size_t size() {
    return _size/sizeof(Scalar);
  }
};


}


#endif // __FLUC_LINALG_WORKSPACE_HH__
