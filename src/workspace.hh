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
template<class Scalar>
class Workspace
{
protected:
  /**
   * Holds the size of the allocated memory.
   */
  size_t size;

  /**
   * Holds a pointer to the allocated workspace memory.
   */
  Scalar *data;


public:
  /**
   * Creates a workspace with the given size.
   */
  Workspace(size_t size=0)
    : size(size), data(0)
  {
    // Preallocate some space:
    if (0 < size)
    {
      this->data = new Scalar[size];
    }
  }

  /**
   * Frees the allocated memory.
   */
  ~Workspace()
  {
    // Free workspace if allocated:
    if (0 != this->data)
    {
      delete this->data;
    }
  }

  /**
   * Ensures that at least @c size elements are allocated for the working-memory.
   */
  void ensure(size_t size)
  {
    // If more space is requested than allocated -> allocate more
    if (this->size < size)
    {
      if (0 != this->data)
      {
        delete this->data;
      }

      this->data = new Scalar[size];
      this->size = size;
    }
  }

  /**
   * Returns a pointer to the allocated working memory.
   */
  Scalar* operator *()
  {
    return this->data;
  }
};


}


#endif // __FLUC_LINALG_WORKSPACE_HH__
