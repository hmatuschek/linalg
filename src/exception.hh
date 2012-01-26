/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */


#ifndef __FLUC_LINALG_EXCEPTION_HH__
#define __FLUC_LINALG_EXCEPTION_HH__

#include <exception>
#include <string>
#include <sstream>


#define LINALG_SHAPE_ASSERT(exp) if (!(exp)) { \
  Linalg::IndexError err; \
  err << "ShapeError at " << __FILE__ \
  << " (" << __LINE__ << "): " \
  << #exp << " failed."; \
  throw err; \
  }

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



class MemoryError: public Exception
{
public:
  MemoryError(const std::string &msg = "")
    : Exception(msg)
  {
    // Pass...
  }

  MemoryError(const MemoryError &other)
    : Exception(other)
  {
    // Pass...
  }

  virtual ~MemoryError() throw ()
  {
  }
};



class LapackError : public Exception
{
public:
  LapackError(const std::string &msg="")
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


class PythonError : public Exception
{
public:
  PythonError(const std::string &msg="")
    : Exception(msg)
  {
    // Pass...
  }

  PythonError(const PythonError &other)
    : Exception(other)
  {
    // Pass...
  }

  ~PythonError() throw ()
  {

  }
};

}


#endif // EXCEPTION_HH
