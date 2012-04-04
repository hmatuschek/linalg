/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_NRM2_HH__
#define __LINALG_BLAS_NRM2_HH__


extern "C" {
double dnrm2_(int *N, const double *x, int *INC);
}


#include "blas/utils.hh"
#include "vector.hh"


namespace Linalg {
namespace Blas {


/**
 * Calculates 2-norm of vector.
 *
 * \f[nrm2(x) = \sqrt{dot(x,x)} \f]
 *
 * @ingroup blas1
 */
inline double nrm2(const Vector<double> &x)
{
  int N   = BLAS_DIMENSION(x);
  int INC = BLAS_INCREMENT(x);
  return dnrm2_(&N, x.ptr(), &INC);
}


}
}
#endif // __LINALG_BLAS_NRM2_HH__
