#ifndef POTRFTEST_HH
#define POTRFTEST_HH

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

public:
  virtual void setUp();

  void testRowMajor();
  void testColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // POTRFTEST_HH
