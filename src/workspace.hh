#ifndef __FLUC_LINALG_WORKSPACE_HH__
#define __FLUC_LINALG_WORKSPACE_HH__

#include <cstdlib>


namespace Fluc {
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
  Workspace(size_t size=0);
  ~Workspace();

  void ensure(size_t size);

  double* operator *();
};


}
}

#endif // __FLUC_LINALG_WORKSPACE_HH__
