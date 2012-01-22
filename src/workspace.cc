#include "workspace.hh"

using namespace Fluc::Linalg;


Workspace::Workspace(size_t size)
  : size(size), data(0)
{
  // Preallocate some space:
  if (0 < size)
  {
    this->data = new double[size];
  }
}


Workspace::~Workspace()
{
  // Free workspace if allocated:
  if (0 != this->data)
  {
    delete this->data;
  }
}


void
Workspace::ensure(size_t size)
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


double *
Workspace::operator *()
{
  return this->data;
}
