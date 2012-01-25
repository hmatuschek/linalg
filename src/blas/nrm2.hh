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

#include "vector.hh"


namespace Linalg {
namespace Blas {


/**
 * Calculates 2-norm of vector.
 *
 * \f[\text{nrm2}(x) = \sqrt{\text{dot}(x)}
 *
 * @ingroup blas1
 */
double nrm2(const Vector<double> &x)
{
  int N   = x.dim();
  int INC = x.stride();

  return dnrm2_(&N, *x, &INC);
}


}
}
#endif // __LINALG_BLAS_NRM2_HH__
