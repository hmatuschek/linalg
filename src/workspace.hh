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
  size_t size;
  double *data;


public:
  Workspace(size_t size=0)
    : size(size), data(0)
  {
    // Preallocate some space:
    if (0 < size)
    {
      this->data = new double[size];
    }
  }

  ~Workspace()
  {
    // Free workspace if allocated:
    if (0 != this->data)
    {
      delete this->data;
    }
  }

  void ensure(size_t size)
  {
    // If more space is requested than allocated -> allocate more
    if (this->size < size)
    {
      if (0 != this->data)
      {
        delete this->data;
      }

      this->data = new double[size];
      this->size = size;
    }
  }

  double* operator *()
  {
    return this->data;
  }
};


}


#endif // __FLUC_LINALG_WORKSPACE_HH__
