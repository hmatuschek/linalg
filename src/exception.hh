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
 * @ingroup error
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
  /**
   * Constructs an exception with given matrix.
   */
  Exception(const std::string &msg = "")
    : std::exception()
  {
    (*this) << msg;
  }

  /**
   * Copy constructor.
   */
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
  /**
   * Constructs a MemoryException with given message.
   */
  MemoryError(const std::string &msg = "")
    : Exception(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  MemoryError(const MemoryError &other)
    : Exception(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
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
  /**
   * Constructs an @c RuntimeError exception from given message.
   */
  RuntimeError(const std::string &msg = "")
    : Exception(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  RuntimeError(const RuntimeError &other)
    : Exception(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
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
  /**
   * Constructs an @c IndexError exception from message.
   */
  IndexError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  IndexError(const IndexError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
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
  /**
   * Constructs a @c ShapeError exception from given message.
   */
  ShapeError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  ShapeError(const ShapeError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
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
  /**
   * Constructs a @c SingularMatrixError exception from message.
   */
  SingularMatrixError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  SingularMatrixError(const SingularMatrixError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
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
  /**
   * Constructs a @c IndefiniteMatrixError exception from message.
   */
  IndefiniteMatrixError(const std::string &msg = "")
    : RuntimeError(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  IndefiniteMatrixError(const IndefiniteMatrixError &other)
    : RuntimeError(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
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
  /**
   * Constructs a @c LapackError exception from given message.
   */
  LapackError(const std::string &msg="")
    : Exception(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  LapackError(const LapackError &other)
    : Exception(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
  ~LapackError() throw ()
  {
  }
};



/**
 * This is the base class of exceptions thrown by the Python interfaces.
 */
class PythonError : public Exception
{
public:
  /**
   * Constructs a @c PythonError from given message.
   */
  PythonError(const std::string &msg="")
    : Exception(msg)
  {
    // Pass...
  }

  /**
   * Copy constructor.
   */
  PythonError(const PythonError &other)
    : Exception(other)
  {
    // Pass...
  }

  /**
   * Destructor.
   */
  ~PythonError() throw ()
  {
    // pass...
  }
};

}


#endif // EXCEPTION_HH
