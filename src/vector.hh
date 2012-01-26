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
public:
  /**
   * A temporary vector owning the data-pointer.
   *
   * If this temporary vector is not assigned to an other Vector instance, the
   * data will be freed.
   */
  typedef std::auto_ptr< Vector<Scalar> > unowned;

  /**
   * Represents a weak reference to the values of the vector.
   */
  class ValueRef
  {
  protected:
    /**
     * Holds a weak reference to the Vector.
     */
    Vector<Scalar> _vector;

  public:
    ValueRef(const Vector &vector)
      : _vector(vector)
    {
      // Pass...
    }

    ValueRef(const ValueRef &other)
      : _vector(other._vector)
    {
      // Pass...
    }

    const ValueRef &operator= (const Scalar &value)
    {
      for (size_t i=0; i<this->_vector.dim(); i++) {
        this->_vector(i) = value;
      }
    }

    const ValueRef &operator= (const Vector<Scalar> &other)
    {
      LINALG_SHAPE_ASSERT(this->_vector.dim() == other.dim());
      for (size_t i=0; i<this->_vector.dim(); i++) {
        this->_vector(i) = other(i);
      }
    }
  };


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
  Vector()
    : DataPtr<Scalar>(), _dimension(0), _offset(0), _increment(0)
  {
    // Pass...
  }

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

  explicit Vector(unowned other)
    : DataPtr<Scalar>(*other, true), _dimension(other->_dimension), _offset(other->_offset),
      _increment(other->_increment)
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
  inline unowned copy()
  {
    Vector out(this->dim());

    for (size_t i=0; i<this->dim(); i++)
      out(i) = (*this)(i);

    return out.takeOwnership();
  }


  inline const ValueRef &vals()
  {
    return ValueRef(*this);
  }


  inline unowned takeOwnership()
  {
    return unowned(new Vector<Scalar>(*this, true));
  }


  inline const Vector<Scalar> operator= (const Vector<Scalar> &other)
  {
    if (this->_owns_data && this->_data != other._data) {
      delete this->_data;
    }

    this->_data = other._data;
    this->_owns_data = false;

    this->_dimension = other._dimension;
    this->_offset = other._offset;
    this->_increment = other._increment;

    return *this;
  }


  inline const Vector<Scalar> operator= (unowned &other)
  {
    if (this->_owns_data && this->_data != other->_data) {
      delete this->_data;
    }

    this->_data = other->_data; other->releaseData();
    this->_owns_data = true;

    this->_dimension = other->_dimension;
    this->_offset = other->_offset;
    this->_increment = other->_increment;

    return *this;
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
};


}

#endif
