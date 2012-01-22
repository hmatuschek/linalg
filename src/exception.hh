#ifndef __FLUC_LINALG_EXCEPTION_HH__
#define __FLUC_LINALG_EXCEPTION_HH__

#include <exception.hh>


namespace Fluc {
namespace Linalg {

class Exception : public Fluc::Exception
{
public:
  Exception(const std::string &msg = "")
    : Fluc::Exception("Linalg: ")
  {
    (*this) << msg;
  }

  ~Exception() throw()
  {}

  Exception(const Exception &other)
    : Fluc::Exception(other)
  {
    // Pass...
  }
};



class IndexError: public Exception
{
public:
  IndexError(const std::string &msg = "")
    : Exception(msg)
  {
    // Pass...
  }

  ~IndexError() throw ()
  {

  }

  IndexError(const IndexError &other)
    : Exception(other)
  {
    // Pass...
  }
};



class LapackError : public Exception
{
public:
  LapackError(const std::string &msg)
    : Exception(msg)
  {
    // Pass...
  }

  LapackError(const LapackError &other)
    : Exception(other)
  {
    // Pass...
  }

  ~LapackError() throw ()
  {

  }
};
}
}

#endif // EXCEPTION_HH
