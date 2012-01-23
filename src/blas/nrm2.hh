#ifndef __LINALG_BLAS_NRM2_HH__
#define __LINALG_BLAS_NRM2_HH__


extern "C" {
void dnrm2_(int *N, double *x, int *INC);
}


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
