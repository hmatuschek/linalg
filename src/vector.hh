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
protected:
  /**
   * Assembles a complete vector from given memory.
   */
  Vector(Scalar *data, size_t offset, size_t dim, size_t incr, bool takes_ownership=true)
    : Array<Scalar>(DataPtr<Scalar>(new DataMngr<Scalar>(data, takes_ownership)),
        offset, std::vector<size_t>(1, dim), std::vector<size_t>(1, incr))
  {
    // Pass...
  }


public:
  /**
   * Assembles a complete vector from given memory.
   */
  Vector(const DataPtr<Scalar> &data, size_t offset, size_t dim, size_t incr)
    : Array<Scalar>(data, offset, std::vector<size_t>(1, dim), std::vector<size_t>(1, incr))
  {
    // Pass...
  }


  /**
   * Constructs an empty vector.
   */
  Vector()
    : Array<Scalar>()
  {
    // Pass...
  }


  Vector(size_t dim)
    : Array<Scalar>()
  {
    Scalar *data = new Scalar[dim];
    *this = Vector<Scalar>(data, 0, dim, 1, true);
  }


  /**
   * Copy constructor, does not copy the memory.
   */
  Vector(const Vector &other)
    : Array<Scalar>(other)
  {
    // Pass...
  }

  /**
   * Copy constructor from array.
   */
  explicit Vector(const Array<Scalar> &other)
    : Array<Scalar>(other)
  {
    LINALG_SHAPE_ASSERT(1 == other.ndim());
  }

  /**
   * Copy constructor, also copies the memory.
   */
  inline Vector<Scalar> copy()
  {
    return Vector<Scalar>(Array<Scalar>::copy());
  }


  /**
   * Assignement of an other vector (weak reference).
   */
  inline Vector<Scalar> &operator= (const Vector<Scalar> &other)
  {
    Array<Scalar>::operator= (other);
    return *this;
  }


  /**
   * Returns the dimension of the vector.
   */
  inline size_t dim() const
  {
    return Array<Scalar>::shape(0);
  }


  /**
   * Returns the offset within memory.
   */
  inline size_t offset() const
  {
    return Array<Scalar>::offset();
  }


  /**
   * Returns the increment within memory.
   */
  inline size_t stride() const
  {
    return Array<Scalar>::strides(0);
  }


  /**
   * Returns a reference to the i-th element.
   */
  inline Scalar &operator ()(size_t i)
  {
    return this->_data[offset() + i*this->stride()];
  }


  /**
   * Returns a const reference to the i-th element.
   */
  inline const Scalar &operator ()(size_t i) const
  {
    return this->_data[offset() + i*this->stride()];
  }


  /**
   * Returns a sub-vector of this vector.
   */
  inline Vector<Scalar> sub(size_t i, size_t n)
  {
    return Vector<Scalar>(*this, offset()+stride()*i, n, stride());
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
