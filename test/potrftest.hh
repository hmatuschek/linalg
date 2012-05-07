#ifndef POTRFTEST_HH
#define POTRFTEST_HH

#include <complex>
#include "unittest.hh"
#include "matrix.hh"
#include "symmatrix.hh"
#include "trimatrix.hh"


class POTRFTest : public UnitTest::TestCase
{
private:
  Linalg::Matrix<double> cA;
  Linalg::Matrix<double> fA;
  Linalg::Matrix<double> cACholU;

  Linalg::Matrix<double> A, Acu, Acl;
  Linalg::Matrix< std::complex<double> > B, Bcu, Bcl;

public:
  virtual void setUp();

  void testDPOTRFRowMajor();
  void testDPOTRFColMajor();

  void testPOTRFRealLower();
  void testPOTRFRealUpper();

  void testPOTRFCmplxLower();
  void testPOTRFCmplxUpper();

public:
  static UnitTest::TestSuite *suite();
};

#endif // POTRFTEST_HH
