#ifndef __LINLAG_CHOLESKYTEST_HH__
#define __LINALG_CHOLESKYTEST_HH__


#include "unittest.hh"
#include "cholesky.hh"


class CholeskyTest : public UnitTest::TestCase
{
protected:
  Linalg::Matrix<double> A, Acu, Acl;
  Linalg::Matrix< std::complex<double> > B, Bcu, Bcl;

public:
  virtual void setUp();

  void testRealLower();
  void testRealUpper();

  void testCmplxLower();
  void testCmplxUpper();

public:
  static UnitTest::TestSuite *suite();
};

#endif // CHOLESKYTEST_HH
