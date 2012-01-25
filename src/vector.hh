#ifndef __LINALG_VECTOR_HH__
#define __LINALG_VECTOR_HH__

#include "memory.hh"
#include "exception.hh"


namespace Linalg {

/**
 * Implements a view to an array (vector).
 *
 * @ingroup core
 */
template <class Scalar>
class Vector : public DataPtr<Scalar>
{
protected:
  /**
   * Holds the dimension (number of elements) of the vector.
   */
  size_t _dimension;

  /**
   * Offset in sizeof(T) within the memory.
   */
  size_t _offset;

  /**
   * Increment between two vector elements (in sizeof(T)), allows to address rows within
   * matrices.
   */
  size_t _increment;


public:
  /**
   * Assembles a complete vector from given memory.
   */
  Vector(Scalar *data, size_t dim, size_t offset, size_t incr, bool takes_ownership)
    : DataPtr<Scalar>(data, takes_ownership), _dimension(dim), _offset(offset), _increment(incr)
  {
    // Pass...
  }

  /**
   * Copy constructor, does not copy the memory.
   */
  explicit Vector(Vector &other, bool take_ownership)
    : DataPtr<Scalar>(other, take_ownership), _dimension(other._dimension), _offset(other._offset),
      _increment(other._increment)
  {
    // Pass...
  }

  /**
   * Const copy constructor.
   */
  Vector(const Vector &other)
    : DataPtr<Scalar>(other), _dimension(other._dimension), _offset(other._offset),
      _increment(other._increment)
  {
    // Pass...
  }

  /**
   * Constructs a new vector with the given dimension.
   */
  Vector(size_t dim)
    : DataPtr<Scalar>(dim), _dimension(dim), _offset(0), _increment(1)
  {
  }


  /**
   * Copy constructor, also copies the memory.
   */
  Vector copy()
  {
    Vector out(this->dim());

    for (size_t i=0; i<this->dim(); i++)
      out(i) = (*this)(i);
  }


  /**
   * Returns the dimension of the vector.
   */
  inline size_t dim() const
  {
    return this->_dimension;
  }


  /**
   * Returns the offset within memory.
   */
  inline size_t offset() const
  {
    return this->_offset;
  }


  /**
   * Returns the increment within memory.
   */
  inline size_t stride() const
  {
    return this->_increment;
  }


  /**
   * Retunrs a pointer to the first element.
   */
  inline Scalar *operator *()
  {
    return this->_data + this->offset();
  }


  /**
   * Returns a const pointer to the first element.
   */
  inline const Scalar *operator *() const
  {
    return this->_data + this->offset();
  }


  /**
   * Returns a reference to the i-th element.
   */
  inline Scalar &operator ()(size_t i)
  {
    return this->_data[this->_offset + i*this->_increment];
  }


  /**
   * Returns a const reference to the i-th element.
   */
  inline const Scalar &operator ()(size_t i) const
  {
    return this->_data[this->_offset + i*this->_increment];
  }

  /**
   * Returns a sub-vector of this vector.
   */
  inline Vector<Scalar> sub(size_t i, size_t n)
  {
    return Vector<Scalar>(this->_data,
                          n, this->offset()+this->_increment*i, this->_increment,
                          false);
  }


  void set(const Vector<Scalar> &other)
  {
    if (! this->dim() == other.dim())
    {
      Linalg::IndexError err;
      err << "Shape mismatch!";
      throw err;
    }

    for (size_t i=0; i<this->dim; i++)
    {
      (*this)(i) = other(i);
    }
  }


  inline void set(size_t i, size_t n, const Vector<Scalar> &other)
  {
    this->sub(i,n).set(other);
  }
};


}

#endif
