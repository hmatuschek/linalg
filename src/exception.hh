/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

/**
 * @defgroup error Error handling and exceptions
 */


#ifndef __FLUC_LINALG_EXCEPTION_HH__
#define __FLUC_LINALG_EXCEPTION_HH__

#include <exception>
#include <string>
#include <sstream>


/**
 * Simple helper macro, that throws an @c IndexError exception if some assertion is not met.
 *
 * @ingroup error.
 */
#define LINALG_SHAPE_ASSERT(exp) if (!(exp)) { \
  Linalg::ShapeError err; \
  err << "ShapeError at " << __FILE__ \
  << " (" << __LINE__ << "): " \
  << #exp << " failed."; \
  throw err; \
  }



namespace Linalg {


/**
 * Base class of all exceptions.
 *
 * This class is derived from @c std::exception, that defines a common interface for all exceptions
 * thrown in C++. However, this class also derives from @c std::ostringstream, that allows to
 * assemble an error message using the "streaming operator" <<.
 *
 * @ingroup error
 */
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



/**
 * Will be thrown if there is an error in the memory-mangement.
 *
 * @ingroup error
 */
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



/**
 * Base class for all runtime-exceptions.
 *
 * @ingroup error
 */
class RuntimeError : public Exception
{
public:
  RuntimeError(const std::string &msg = "")
    : Exception(msg)
  {
    // Pass...
  }

  RuntimeError(const RuntimeError &other)
    : Exception(other)
  {
    // Pass...
  }

  virtual ~RuntimeError() throw ()
  {
  }
};



/**
 * Will be thrown, if an index is out-of-range or some shape-assertion is not met.
 *
 * @ingroup error
 */
class IndexError: public RuntimeError
{
public:
  IndexError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  IndexError(const IndexError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  virtual ~IndexError() throw ()
  {
  }
};



/**
 * Will be thrown if some shape-assertion is not met.
 *
 * @ingroup error
 */
class ShapeError: public RuntimeError
{
public:
  ShapeError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  ShapeError(const ShapeError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  virtual ~ShapeError() throw ()
  {
  }
};



/**
 * Will be thrown if a matrix is singular.
 *
 * @ingroup error
 */
class SingularMatrixError : public RuntimeError
{
public:
  SingularMatrixError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  SingularMatrixError(const SingularMatrixError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  virtual ~SingularMatrixError() throw ()
  {
  }
};



/**
 * Will be thrown, if a matrix is not definite.
 *
 * @ingroup error
 */
class IndefiniteMatrixError : public RuntimeError
{
public:
  IndefiniteMatrixError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  IndefiniteMatrixError(const IndefiniteMatrixError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  virtual ~IndefiniteMatrixError() throw ()
  {
  }
};



/**
 * Will be thrown if an argument to a LAPACK function is invalid.
 *
 * @ingroup error
 */
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
