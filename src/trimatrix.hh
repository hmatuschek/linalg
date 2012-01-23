#ifndef __LINALG_TRIMATRIX_HH__
#define __LINALG_TRIMATRIX_HH__

#include "matrix.hh"

namespace Linalg {

/**
 * Implements a view on @c Matrix as a unit or non-unit upper- or lower-triangular matrix.
 */
template <class Scalar>
class TriMatrix : public Matrix<Scalar>
{
protected:
  bool _is_upper;
  bool _is_unit_triangular;


public:
  TriMatrix(Matrix<Scalar> &matrix, bool upper, bool unit)
    : Matrix<Scalar>(matrix.weak()), _is_upper(upper), _is_unit_triangular(unit)
  {
    if (matrix.transposed())
      this->_is_upper = ! this->_is_upper;
  }


  TriMatrix(TriMatrix<Scalar> &other)
    : Matrix<Scalar>(other.weak()),
      _is_upper(other._is_upper), _is_unit_triangular(other._is_unit_triangular)
  {
    // Pass...
  }


  inline bool isUpper() const
  {
    if (this->_transposed)
      return ! this->_is_upper;
    return this->_is_upper;
  }


  inline bool isUnit() const
  {
    return this->_is_unit_triangular;
  }

};

}


#endif // __LINALG_TRIMATRIX_HH__
