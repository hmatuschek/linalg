#ifndef __LINALG_VECTOR_HH__
#define __LINALG_VECTOR_HH__

#include "array.hh"
#include "exception.hh"


namespace Linalg {

/**
 * Implements a view to an array (vector).
 *
 * @ingroup matrix
 */
template <class Scalar>
class Vector
    : public Array<Scalar>
{
public:
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
    /**
     * Constructs a @c ValueRef (a reference to the values of a vector) for the given vector.
     */
    ValueRef(const Vector &vector)
      : _vector(vector)
    {
      // Pass...
    }

    /**
     * Copy constructor.
     */
    ValueRef(const ValueRef &other)
      : _vector(other._vector)
    {
      // Pass...
    }

    /**
     * Assignment operation for a scalar.
     */
    const ValueRef &operator= (const Scalar &value)
    {
      for (size_t i=0; i<this->_vector.dim(); i++) {
        this->_vector(i) = value;
      }
    }


    /**
     * Assignment operation for an other vector.
     */
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
   * Assembles a complete vector from given memory.
   */
  Vector(Scalar *data, size_t offset, size_t dim, size_t incr, bool takes_ownership=true)
    : Array<Scalar>(data, _offset(offset), _shape(1, dim), _strides(1, incr), takes_ownership)
  {
    // Pass...
  }

  /**
   * Assembles a complete vector from given memory.
   */
  Vector(DataPtr<Scalar> *data, size_t offset, size_t dim, size_t incr)
    : Array<Scalar>(data, _offset(offset), _shape(1, dim), _strides(1, incr))
  {
    // Pass...
  }


public:
  /**
   * Constructs an empty vector.
   */
  Vector()
    : Array()
  {
    // Pass...
  }

  /**
   * Copy constructor, does not copy the memory.
   */
  Vector(const Vector &other)
    : Array(other)
  {
    // Pass...
  }


  /**
   * Copy constructor, also copies the memory.
   */
  inline Vector<Scalar> copy()
  {
    return Vector<Scalar>(Array<Scalar>::copy());
  }


  /**
   * Returns a reference to the values of this vector.
   */
  inline ValueRef vals()
  {
    return ValueRef(*this);
  }


  /**
   * Assignement of an other vector (weak reference).
   */
  inline Vector<Scalar> &operator= (const Vector<Scalar> &other)
  {
    Array<Scalar>::operator =(other);
    return *this;
  }


  /**
   * Takes the ownership of the unowned vector.
   */
  inline const Vector<Scalar> &operator= (unowned other)
  {
    Array<Scalar>::operator =(other);

    return *this;
  }


  /**
   * Returns the dimension of the vector.
   */
  inline size_t dim() const
  {
    return _shape[0];
  }


  /**
   * Returns the offset within memory.
   */
  inline size_t offset() const
  {
    return _offset;
  }


  /**
   * Returns the increment within memory.
   */
  inline size_t stride() const
  {
    return _strides[0];
  }


  /**
   * Retunrs a pointer to the first element.
   */
  inline Scalar *operator *()
  {
    return this->_data->ptr() + this->offset();
  }


  /**
   * Returns a const pointer to the first element.
   */
  inline const Scalar *operator *() const
  {
    return this->_data->ptr() + this->offset();
  }


  /**
   * Returns a reference to the i-th element.
   */
  inline Scalar &operator ()(size_t i)
  {
    return this->_data->ptr()[_offset + i*this->stride()];
  }


  /**
   * Returns a const reference to the i-th element.
   */
  inline const Scalar &operator ()(size_t i) const
  {
    return this->_data->ptr()[_offset + i*this->stride()];
  }


  /**
   * Returns a sub-vector of this vector.
   */
  inline Vector<Scalar> sub(size_t i, size_t n)
  {
    return Vector<Scalar>(_data, offset()+stride()*i, n, stride());
  }


public:
  /**
   * Allocates an uninitialized vector of given size.
   */
  static Vector<Scalar> empty(size_t dim)
  {
    Scalar *data = new Scalar[dim];
    return Vector(data, 0, dim, 1, true);
  }


  static Vector<Scalar> zero(size_t dim)
  {
    Vector<Scalar> vec = empty(dim);
    for (size_t i=0; i<dim; i++)
      vec(i) = 0;
    return vec;
  }
};


}

#endif
