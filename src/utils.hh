#ifndef __LINALG_UTILS_HH__
#define __LINALG_UTILS_HH__

#include <complex>


namespace Linalg {

/**
 * Utility function to unify the Cholesky decomposition codes. This (inline) function calculates
 * a*b if a,b are real and a*b^* if a,b are complex.
 */
inline double
__prod_abcc(const double &a, const double &b)
{
  return a*b;
}

/**
 * Utility function to unify the Cholesky decomposition codes. This (inline) function calculates
 * a*b if a,b are real and a*b^* if a,b are complex.
 */
inline std::complex<double>
__prod_abcc(const std::complex<double> &a, const std::complex<double> &b)
{
  return std::complex<double>(a.real()*b.real() + a.imag()*b.imag(),
                              a.imag()*b.real() - a.real()*b.imag());
}

/**
 * Utility function to calculate either a*a if a is real or a*a^* if a is complex.
 */
inline double
__prod_aacc(const double &a) {
  return a*a;
}

/**
 * Utility function to calculate either a*a if a is real or a*a^* if a is complex.
 */
inline double
__prod_aacc(const std::complex<double> &a) {
  return a.real()*a.real() + a.imag()*a.imag();
}


/**
 * Utility function to test if a value is positive and real.
 */
inline bool __is_real_pos(const double &a)
{
  return a > 0;
}

/**
 * Utility function to test if a value is positive and real.
 */
inline bool __is_real_pos(const std::complex<double> &a) {
  return a.imag() == 0 && a.real() > 0;
}


/**
 * Utility function to get the real part of a complex or real value.
 */
inline double __get_real(const double &a) {
  return a;
}

/**
 * Utility function to get the real part of a complex or real value.
 */
inline double __get_real(const std::complex<double> &a) {
  return a.real();
}

/**
 * Utility function to get the imag part of a complex or real value.
 */
inline double __get_imag(const double &a) {
  return 0.0;
}

/**
 * Utility function to get the imag part of a complex or real value.
 */
inline double __get_imag(const std::complex<double> &a) {
  return a.imag();
}


}


#endif // __LINALG_UTILS_HH__
