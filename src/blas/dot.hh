#ifndef __LINALG_BLAS_DOT_HH__
#define __LINALG_BLAS_DOT_HH__

#include "vector.h"

extern "C" {
ddot_(const int *N, const double *DX, const int *INCX, double *DY, const int *INCY);
}


namespace Linalg {
namespace Blas {


/**
 * Calculates simple dot product:
 *
 * \f[dot(x,y) = x^T\cdot y\f]
 *
 * @ingroup blas1
 */
double dot(const Vector<double> &x, const Vector<double> &y) {
  if (x.dim() != y.dim()) {
    IndexError err; err << "Dimension mismatch."; throw err;
  }

  int N = x.dim();
  int INCX = x.stride();
  int INCY = y.stride();

  return ddot_(&N, *x, &INCX, *y, &INCY);
}


}
}
#endif // __LINALG_BLAS_DOT_HH__
