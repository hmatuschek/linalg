#ifndef __FLUC_LINALG_EXCEPTION_HH__
#define __FLUC_LINALG_EXCEPTION_HH__

#include <exception>
#include <string>
#include <sstream>


namespace Linalg {

class Exception : public std::exception, public std::ostringstream
{
public:
  Exception(const std::string &msg = "")
    : std::exception()
  {
    (*this) << msg;
  }

  Exception(const Exception &other)
    : std::exception()
  {
    (*this) << other.what();
  }

  /**
   * Is needed by the std::exception interface. Does nothing.
   */
  virtual ~Exception() throw ()
  {
    // Pass...
  }

  /**
   * Returns the message of the exception. This method is needed by the std::exception interface.
   */
  virtual const char *what() const throw()
  {
    return this->str().c_str();
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

  IndexError(const IndexError &other)
    : Exception(other)
  {
    // Pass...
  }

  virtual ~IndexError() throw ()
  {
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


#endif // EXCEPTION_HH
