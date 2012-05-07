#ifndef GEQRFTEST_HH
#define GEQRFTEST_HH

#include "unittest.hh"
#include <matrix.hh>
#include <vector.hh>


class GEQRFTest : public UnitTest::TestCase
{
private:
  Linalg::Matrix<double> A;
  Linalg::Vector<double> tau;

public:
  virtual void setUp();

  void testRealNN();
  void testRealMN();

public:
  static UnitTest::TestSuite *suite();
};

#endif // QRTEST_HH
