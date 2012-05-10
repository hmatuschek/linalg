#ifndef __LINALG_VECTOR_OPERATORS_HH__
#define __LINALG_VECTOR_OPERATORS_HH__

#include "blas/scal.hh"
#include "blas/nrm2.hh"


namespace Linalg {


/**
 * Implements the in-place vector-scaleing.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Vector<Scalar> & operator *= (Vector<Scalar> &V, const Scalar &alpha) {
  Blas::scal(alpha, V);
  return V;
}




/**
 * Implements vector-scaleing \f$y' = \frac{1}{\alpha}y$.
 *
 * @ingroup operators
 */
template <class Scalar>
inline Vector<Scalar> & operator /= (Vector<Scalar> &V, const Scalar &alpha) {
  Blas::scal(Scalar(1)/alpha, V);
  return V;
}


/**
 * Streaming operator for vectors.
 *
 * @ingroup operators
 */
template <class Scalar>
std::ostream &operator<< (std::ostream &stream, Vector<Scalar> &vector)
{
  stream << "[";
  for (size_t i=0; i<vector.dim()-1; i++) {
    stream << vector(i) << ", ";
  }
  stream << vector(vector.dim()-1) << "]";

  return stream;
}

}



namespace std {

/**
 * Extends the @c std::abs2 function to hanlde vectors as \f$abs2(v) = v^Tv\f$
 *
 * @ingroup operators
 */
inline double abs2(const Linalg::Vector<double> &v)
{
  return Linalg::Blas::nrm2sq(v);
}


/**
 * Extends the @c std::abs function to hanlde vectors as \f$abs(v) = \sqrt(v^Tv)\f$
 *
 * @ingroup operators
 */
inline double abs(const Linalg::Vector<double> &v) {
  return Linalg::Blas::nrm2(v);
}

}


#endif // __LINALG_VECTOR_OPERATORS_HH__
